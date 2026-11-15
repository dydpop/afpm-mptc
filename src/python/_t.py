"""
表贴式 PMSM + MPTC / FOC 数值仿真（纯 numpy，不依赖 Simulink）

用途：
1. 在搭 Simulink 之前，先验证 MPTC 算法逻辑对不对
2. 参数扫描时快速出数据（比 Simulink 快 100 倍）
3. 给组员搭模型提供"标准答案"对照

约定：
- αβ 静止坐标系，Clark 为等幅值变换
- 电机为表贴式（Ld = Lq = L），无磁阻转矩
- 逆变器两电平，8 种开关状态
- 控制器采样 Ts = 100 μs，数值积分步长 Ts/10
"""

import numpy as np

# ============ 电机与系统参数（文献 Siami 2013 轴向磁通样机） ============
P = dict(
    Rs=0.2,          # 定子电阻 Ω
    L=8.5e-3,        # d/q 电感 H（表贴式 Ld=Lq）
    psi=0.175,       # 永磁磁链 Wb
    p=4,             # 极对数
    J=0.003,         # 转动惯量 kg·m²
    F=0.008,         # 阻尼 N·m·s
    Udc=250.0,       # 母线电压 V
    Ts=50e-6,        # 控制器采样周期（20 kHz）
    Nsub=10,         # 每个 Ts 内的积分子步数
)

# 额定值（用于代价函数归一化）
T_RATED = 11.9       # N·m，按 0.5 HP / 300 r/min 折算：373 W / 31.4 rad·s⁻¹
W_RATED = 300 * 2*np.pi/60   # rad/s ≈ 31.4


def voltage_vectors(Udc):
    """返回 8 个开关状态对应的 (u_alpha, u_beta) 以及开关状态本身"""
    states = np.array([
        [1,0,0],[1,1,0],[0,1,0],[0,1,1],[0,0,1],[1,0,1],  # 6 个非零
        [0,0,0],[1,1,1],                                   # 2 个零矢量
    ], dtype=float)
    sa, sb, sc = states[:,0], states[:,1], states[:,2]
    ua = (2/3) * Udc * (sa - 0.5*sb - 0.5*sc)
    ub = (2/3) * Udc * (np.sqrt(3)/2) * (sb - sc)
    return ua, ub, states


def deriv(x, ua, ub, prm, T_load):
    """状态导数。x = [i_a, i_b, w_m, theta_m]"""
    ia, ib, wm, thm = x
    Rs, L, psi, p, J, F = prm['Rs'], prm['L'], prm['psi'], prm['p'], prm['J'], prm['F']
    the = p * thm
    # 反电动势
    # psi_a = L*ia + psi*cos(th);  u = R*i + dpsi_a/dt
    # => L*dia/dt = u - R*ia + psi*p*wm*sin(th)
    ea =  psi * p * wm * np.sin(the)
    eb = -psi * p * wm * np.cos(the)
    dia = (ua - Rs*ia + ea) / L
    dib = (ub - Rs*ib + eb) / L
    # 磁链与转矩
    pa = L*ia + psi*np.cos(the)
    pb = L*ib + psi*np.sin(the)
    Te = 1.5 * p * (pa*ib - pb*ia)
    dwm = (Te - T_load - F*wm) / J
    return np.array([dia, dib, dwm, wm]), Te


def flux_torque(ia, ib, thm, prm):
    """由电流和转子位置算磁链幅值与转矩"""
    L, psi, p = prm['L'], prm['psi'], prm['p']
    the = p * thm
    pa = L*ia + psi*np.cos(the)
    pb = L*ib + psi*np.sin(the)
    Te = 1.5 * p * (pa*ib - pb*ia)
    return np.hypot(pa, pb), Te, pa, pb


def predict(ia, ib, ua, ub, wm, thm, prm, h, nsub=10):
    """预测下一拍状态。内部用 nsub 个子步，减小大 Ts 下的离散化误差"""
    Rs, L, psi, p = prm['Rs'], prm['L'], prm['psi'], prm['p']
    hs = h / nsub
    for _ in range(nsub):
        the = p * thm
        ea =  psi * p * wm * np.sin(the)
        eb = -psi * p * wm * np.cos(the)
        ia = ia + hs*(ua - Rs*ia + ea)/L
        ib = ib + hs*(ub - Rs*ib + eb)/L
        thm = thm + hs*wm
    return ia, ib, thm


def run(mode='mptc', mismatch=None, w_ref=300.0, T_load=8.0, t_end=1.2,
        t_step=0.3, lam_flux=0.2, lam_sw=0.05, seed=0,
        record=True):
    """
    主仿真

    mode      : 'mptc' | 'foc'
    mismatch  : dict，控制器认为的参数与实际参数的偏差，如 {'Rs':1.2,'L':0.9,'psi':0.9}
                （电机实际参数始终是 P，改的是控制器内部用的值）
    w_ref     : 转速给定 r/min
    T_load    : 负载转矩 N·m
    t_step    : 转矩给定阶跃时刻 s（之前为空载 0）
    """
    prm = dict(P)
    ctrl_prm = dict(P)
    if mismatch:
        for k, v in mismatch.items():
            ctrl_prm[k] = P[k] * v

    Udc = prm['Udc']
    Ts = prm['Ts']
    Nsub = prm['Nsub']
    h = Ts / Nsub
    N = int(t_end / Ts)

    ua_v, ub_v, states = voltage_vectors(Udc)

    # 状态初始化
    x = np.array([0.0, 0.0, 0.0, 0.0])
    w_ref_rad = w_ref * 2*np.pi/60

    # 控制器状态
    prev_vec = 0
    prev_vec_applied = 6          # 上一拍实际施加的矢量索引（初始给零矢量）
    int_w = 0.0                    # 转速环积分
    int_d = 0.0; int_q = 0.0       # 电流环积分

    # 记录
    rec = dict(t=[], Te=[], wm=[], id=[], iq=[], ia=[], ib=[], vec=[], Te_ref=[])

    sw_count = 0
    prev_state_bits = states[6].copy()   # 初始为 000

    for k in range(N):
        t = k * Ts
        ia, ib, wm, thm = x
        the = prm['p'] * thm
        T_now = T_load if t >= t_step else 0.0
        # 负载在阶跃时刻突加，用于测转矩响应

        # ---------- 控制器 ----------
        if mode == 'mptc':
            # 0) 转速外环：PI 产生转矩给定（含抗积分饱和）
            w_err = w_ref_rad - wm
            Kp_w, Ki_w = 8.0, 40.0
            T_max = 1.5 * T_RATED
            T_raw = Kp_w*w_err + int_w
            T_ref = float(np.clip(T_raw, -T_max, T_max))
            if T_raw == T_ref:                      # 仅在未饱和时积分
                int_w = float(np.clip(int_w + Ki_w*w_err*Ts, -T_max, T_max))
            # 1) 延时补偿：用上一拍实际施加的矢量，从当前状态推进一拍
            ua_p, ub_p = ua_v[prev_vec_applied], ub_v[prev_vec_applied]
            ia1, ib1, thm1 = predict(ia, ib, ua_p, ub_p, wm, thm, ctrl_prm, Ts)
            # 2) 从 k+1 状态出发，对 8 个矢量预测 k+2
            psi_s_ref = ctrl_prm['psi']
            costs = np.zeros(8)
            Te_pred = np.zeros(8)
            for i in range(8):
                ia2, ib2, thm2 = predict(ia1, ib1, ua_v[i], ub_v[i], wm, thm1, ctrl_prm, Ts)
                _, Te2, _, _ = flux_torque(ia2, ib2, thm2, ctrl_prm)
                ps2, _, _, _ = flux_torque(ia2, ib2, thm2, ctrl_prm)
                Te_pred[i] = Te2
                n_sw = int(np.sum(states[i] != prev_state_bits))
                costs[i] = (abs(T_ref - Te2)/T_RATED
                            + lam_flux * abs(psi_s_ref - ps2)/ctrl_prm['psi']
                            + lam_sw * n_sw)
            sel = int(np.argmin(costs))
            ua_cmd, ub_cmd = ua_v[sel], ub_v[sel]
            prev_vec_applied = sel
            sw_count += int(np.sum(states[sel] != prev_state_bits))
            prev_state_bits = states[sel].copy()
            ua_app, ub_app = ua_cmd, ub_cmd

        else:  # FOC
            # 转速 PI → q 轴电流给定
            w_err = w_ref_rad - wm
            Kp_w, Ki_w = 0.5, 20.0
            iq_ref = Kp_w*w_err + int_w
            iq_ref = np.clip(iq_ref, -20.0, 20.0)
            int_w = np.clip(int_w + Ki_w*w_err*Ts, -20.0, 20.0)
            id_ref = 0.0
            # Park
            id_ = ia*np.cos(the) + ib*np.sin(the)
            iq_ = -ia*np.sin(the) + ib*np.cos(the)
            Kp_i, Ki_i = 8.0, 800.0
            ud = Kp_i*(id_ref-id_) + int_d
            uq = Kp_i*(iq_ref-iq_) + int_q
            int_d = np.clip(int_d + Ki_i*(id_ref-id_)*Ts, -Udc, Udc)
            int_q = np.clip(int_q + Ki_i*(iq_ref-iq_)*Ts, -Udc, Udc)
            # 反 Park
            ua_app = np.cos(the)*ud - np.sin(the)*uq
            ub_app = np.sin(the)*ud + np.cos(the)*uq
            # 限幅到逆变器可输出范围
            mag = np.hypot(ua_app, ub_app)
            lim = Udc/np.sqrt(3)
            if mag > lim:
                ua_app *= lim/mag; ub_app *= lim/mag
            # 理想调制，开关次数按 6 次/Ts 计（对应 10 kHz 载波）
            sw_count += 6

        # ---------- 系统积分 ----------
        for _ in range(Nsub):
            dx, Te = deriv(x, ua_app, ub_app, prm, T_now)
            x = x + h * dx

        if record:
            ia, ib, wm, thm = x
            _, Te, _, _ = flux_torque(ia, ib, thm, prm)
            the = prm['p']*thm
            id_ = ia*np.cos(the) + ib*np.sin(the)
            iq_ = -ia*np.sin(the) + ib*np.cos(the)
            rec['t'].append(t); rec['Te'].append(Te); rec['wm'].append(wm)
            rec['id'].append(id_); rec['iq'].append(iq_)
            rec['ia'].append(ia); rec['ib'].append(ib)
            rec['vec'].append(prev_vec_applied if mode=='mptc' else -1)
            rec['Te_ref'].append(T_ref if mode=='mptc' else np.nan)

    out = {k: np.array(v) for k, v in rec.items()}
    out['f_sw'] = sw_count / (6 * t_end)     # 等效平均开关频率
    return out


def metrics(out, t_step=0.3, T_rated=T_RATED, Ts=100e-6):
    """四项指标"""
    t = out['t']; Te = out['Te']
    m = {}
    # 阶跃响应时间：进入 ±5% 带
    i0 = np.searchsorted(t, t_step)
    band = 0.05 * T_rated
    idx = np.where(np.abs(Te[i0:] - T_rated) <= band)[0]
    m['t_rise'] = t[i0 + idx[0]] - t_step if len(idx) else np.nan
    # 稳态脉动峰峰值（最后 0.1 s）
    i_ss = np.searchsorted(t, t[-1] - 0.25)
    m['ripple'] = Te[i_ss:].max() - Te[i_ss:].min()
    # 定子电流 THD：取稳态相电流做 FFT，基波为最大峰
    sig = out['ia'][i_ss:] - out['ia'][i_ss:].mean()
    if len(sig) > 128:
        S = np.abs(np.fft.rfft(sig))
        k0 = int(np.argmax(S[1:])) + 1          # 基波分量（跳过直流）
        if S[k0] > 0:
            harm = np.sqrt(np.sum(S[1:]**2) - S[k0]**2)
            m['thd'] = harm / S[k0] * 100
        else:
            m['thd'] = np.nan
    else:
        m['thd'] = np.nan
    m['f_sw'] = out['f_sw']
    return m


if __name__ == '__main__':
    print("=== 单工况测试：MPTC，无参数偏差 ===")
    o = run(mode='mptc', t_end=0.5)
    mm = metrics(o)
    print(f"  末速 {o['wm'][-1]*60/(2*np.pi):.1f} r/min")
    print(f"  稳态转矩 {o['Te'][-5000:].mean():.2f} N·m（给定 {T_RATED}）")
    for k, v in mm.items():
        print(f"  {k}: {v:.4g}" if isinstance(v, float) else f"  {k}: {v}")
