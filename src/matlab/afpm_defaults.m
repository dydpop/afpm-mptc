function P = afpm_defaults()
% 研究仿真默认输入；并非样机实测参数集。
P.Rs = 0.2;          % 定子电阻，ohm
P.L = 0.0085;        % 表贴式等效电感 Ld=Lq，H
P.psi = 0.175;       % 永磁磁链，Wb
P.p = 4;             % 极对数
P.J = 0.003;         % 研究假设：惯量，kg*m^2，非实测
P.F = 0.008;         % 研究假设：黏性阻尼，N*m*s，非实测
P.Udc = 250;         % 直流母线电压，V
P.Ts = 50e-6;        % 控制周期，s
P.Nsub = 20;         % 每控制周期的真实plant RK4子步数
P.w_ref = 300;       % 机械转速给定，r/min
P.T_load = 8;        % 阶跃负载，N*m；允许为零
P.t_step = 0.3;      % 负载阶跃时刻，s
P.t_end = 2.0;       % 仿真终止时刻，s；速度环慢极点约4 rad/s，1.2s末仍在拖尾
P.T_rated = 11.9;    % 研究采用的额定转矩基准，N*m
P.lam_flux = 1.0;    % 标称候选0.05/0.5/1/2试验选定；非全局最优，扫描前冻结
P.lam_sw = 0.02;     % 研究假设：桥臂翻转代价权重，非实测
P.Kp_w = 2;          % 研究整定：速度PI比例增益，N*m/(rad/s)
P.Ki_w = 8;          % 研究整定：速度PI积分增益，N*m/rad
P.f_sw_foc = 6500;   % 载波设定频率，Hz；不代表实际统计频率
P.delay_comp = true; % 仅控制MPTC预测补偿；物理一拍延时始终保留
P.log_fine = true;   % 默认保留子步数据，用于PWM频谱与指标分析
P.predictor = 'euler'; % MPTC预测离散方式：'euler'或'trapezoid'（显式Heun）
P.foc_gain_mode = 'estimated'; % 'estimated'按估计参数整定（默认）；'nominal'固定标称增益
P.cost_mode = 'weighted';  % MPTC代价函数：'weighted'加权求和（默认）或'lexicographic'分层筛选
P.lex_torque_tol_pct = 10; % 分层第一层转矩容差，%T_rated
P.lex_flux_tol_pct = 10;   % 分层第二层磁链容差，%ψf
P.predict_form = 'current'; % MPTC预测方程：'current'电流微分（默认）或'flux'磁链积分
P.direct_Tref = [];  % 非空时旁路速度环：struct('t',阶跃时刻s,'from',,'to',) N*m
end
