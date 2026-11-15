# 模型目录

一人一文件，不共用：

| 文件 | 负责人 | 说明 |
|---|---|---|
| `baseline.slx` | Haotian Feng | 基准被控对象，参数替换与坐标约定最终版 |
| `mptc.slx` | Wenkai Liu | MPTC 控制器，含枚举、代价函数、延时补偿 |
| `foc.slx` | Bokun Tian | PI 矢量控制对照 |
| `main.slx` | — | 联调用，只做子系统引用，不改他人模块 |

模型定版后不再改动，改动须在群里说明并更新本表。
