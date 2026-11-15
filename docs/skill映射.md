# 技能映射（AFPM 表贴式 PMSM / MPTC vs PI 矢量控制）

读取范围：`~/.workbuddy/skills/` 与 `~/.codex/skills/` 下与科研相关的 SKILL.md，共 21 个（另 4 个目录无 SKILL.md，已注明）。所有描述均摘自实际 SKILL.md，未读到的写"未找到"。

## 1. 阶段 → 技能 → 检查项 → 触发时机

| 阶段 / 横向能力 | 推荐技能 | 该技能提供的具体检查项或流程 | 触发时机 |
|---|---|---|---|
| P2 直接转矩阶跃 | `results-analysis`（只读审计模式） | 先"锁定比较问题"：哪个方法对哪个基线、主指标是什么、重复单元是什么；未定义"重复单元"则不得下"胜出"或显著性结论 | 阶跃前后转矩超调/恢复时间要下结论前 |
| P2 / P4 阶跃与额定负载复核 | `nature-statistics` + `research-vector-figures`(矢量口径) | 分离"测了什么 / 分析单元 / 推断什么"三问；同一指标在正文与表中出现两种精度、或重叠误差棒被说成"更优"时须排查 | 同一指标多处引用、或需要写"优于"时 |
| P3 离散方式对照 | `results-analysis` | 报告完整统计：mean±std、区间、run/seed 数、效应量、多重比较处理；缺 sample/seed 时明确写 blocker，不得编造 | 离散方式/开关频率口径对比成表时 |
| P5 升速 / P6 弱磁 | `nature-figure`（作图）+ `results-analysis`（证据） | 每张主图必须有 purpose、变量、误差棒含义、caption 要求、解读清单；图回答不了"改变了什么判断"即为装饰图 | 转速-指标曲线成图前先写图契约 |
| P7 Simulink 闭环 | 未找到 | 未找到任何针对 Simulink/MATLAB 模型调试或闭环验证的技能；本阶段只能靠 `6verity` 式数值一致性检查人工执行 | — |
| P8 结题报告 | `results-report` → `paper-self-review` → `spec-driven-docx` | `results-report` 硬门槛：缺 analysis-report / stats-appendix / figure-catalog / provenance 时只写 blocker，不写完整报告；`spec-driven-docx` 按"稿件+要求+模板"出合规 docx 与 PDF | 分析包齐备后再成文；排版阶段 |
| P8 答辩 PPT | `nature-paper2ppt` | 将论文/MD 转为演示结构的流程（详见其 SKILL.md） | 报告定稿后 |
| 横向：数值自查 | `physics-lab-report`（方法借用） | 独立恒等式互检；"数据与理论差得大时先怀疑测量/读数，不改原理、不改公式"，不许悄悄乘系数"修正"；偏差 >0.5% 记为待查项 | 每个 case 出数后立即 |
| 横向：结果审查 | `paper-self-review` + `nature-reviewer` | Claim Audit 输出格式：Claim / Verdict(keep\|weaken\|revise\|remove) / Evidence used / Missing evidence / Overclaim risk / Suggested wording；`nature-reviewer` 的 forensic consistency audit 审算术、指标边界、聚合层级、prose-table 一致性、provenance | 每条结论定稿前 |
| 横向：出图 | `matplotlib`（底座）/ `nature-figure`（规范与 QA）/ `4drawio` 或 `research-vector-figures`（非数据图） | `matplotlib`：OO 接口、constrained_layout、300 dpi、避免 jet、色盲友好；`nature-figure`：5 点契约 + 多面板对齐 1.5 pt 容差 + 碰撞审计脚本 | 数据图 / 系统框图分别走不同技能 |
| 横向：写作边界 | `nature-polishing` | 段落结构问题不能靠"发明内容"修补，只能标出；按 claim/evidence/boundary 优先级处理 | 措辞超出证据范围时 |
| 横向：格式验收 | `6verity`（codex） | 数值必须来自结果记录/结果 JSON，关键指标不得与记录冲突；发现冲突不得自行发明新结果；图表引用、占位符、内部路径泄露检查 | 报告导出前最后一关 |

## 2. 只能帮到报告与图，不能替代数值正确性

- `nature-figure` / `matplotlib` / `4drawio` / `research-vector-figures` 只管"画得对、可审、可读"：它们明确"不重画或不重新解释实测数据"（`research-vector-figures` 原文），不做仿真验证。
- `native-figure` 的 data-integrity gate 要求"不得为让图画得快而删数据"，但它**不会**替你判断转矩脉动口径对不对。
- `paper-self-review` 原文声明：其 claim audit "不替代 `results-analysis` 的统计，也不替代引用核验"，只检查"结论是否正确使用了现有证据"。
- `nature-statistics` 原文："是报告与审查技能，不是统计师重算原始数据的替代品"，缺信息写 `AUTHOR_INPUT_NEEDED`。
- `results-analysis` 的强项是"不编造统计"，但它只在**给定产物**上工作；CSV 本身算错，它只会如实呈现错误。
- `spec-driven-docx` 只核对格式是否符合模板，五道闸门只管"来源/证据/严重度/反查/版本"，不管物理。
- 结论：**数值正确性只能来自 MATLAB/Simulink 侧的独立校验（能量/功率平衡、d-q 方程回代、开关频率计数一致性）**，上述技能全部在此之后再介入。

## 3. 本课题自查清单（12 条）

- [ ] 每条结论都指到具体 CSV 行 / 图号 / 日志路径（`paper-self-review`：Main claims supported by named evidence）。
- [ ] 每条强结论标 Verdict，若为 weaken/remove 必须写出"允许措辞"与"禁止更强措辞"（`results-analysis` claim candidate 模板）。
- [ ] 转矩脉动、THD、等效开关频率三指标各自声明：定义式、统计窗口（稳态段起止）、方向（越小越好）。
- [ ] 同一指标在摘要、正文、表格、图注中只保留一种精度与一种单位（`nature-statistics`：one metric at two precisions）。
- [ ] 参数失配扫描的每个点都注明：改动的是 Rs / L / ψf 中的哪一个、偏离百分比、是否同时联动。
- [ ] 任何"MPTC 优于 PI"的说法必须有 PI 基线同工况数据（`rigor-reviewer` D1：Improvement 类主张需要基线对比）。
- [ ] 若给误差棒/区间，须说明是哪个重复单元（同一 case 多次？不同扫描点？），不得把不同失配点当独立样本做显著性。
- [ ] 闭环未完成时必须明写边界："本文未实现完整 Simulink 闭环，仅数值仿真 + 电机方程校验模型"（`nature-reviewer`：区分 supported / weak / not assessable）。
- [ ] 每张图的 caption 说明误差棒含义、统计窗口、失配条件，正文不得只说"如图 X 所示"（`6verity` 列为问题项）。
- [ ] 图内每个字形 ≥ 5 pt，多面板对齐按 1.5 pt 容差审计，交付前跑碰撞审计（`nature-figure` 契约与 QA）。
- [ ] 派生量（如等效开关频率）与原始数据用独立路径重算一次，偏差 >0.5% 记待查（`physics-lab-report` 步 8）。
- [ ] 报告中的关键数值从结果 JSON/CSV 抽取核对，冲突时回改代码或正文，不得改数（`6verity` Step 5）。

## 4. 不建议使用的技能及原因

| 技能 | 原因 |
|---|---|
| `physics-lab-report` / `physics-lab-report-prep` | 面向大学物理实验课（课本照抄、手写抄写稿 A/B/C 分级、五栏预习）。本课题只**借用**其"独立互检、不改公式、偏差 >0.5% 待查"的方法论，不套用其模板与流程。 |
| `nature-experiment-log` | 面向湿实验日志（马弗炉、管式炉、电化学工作站、样品批次 ID），与纯数值仿真不匹配。 |
| `nature-literature-pipeline` / `research-ideation` / `hypothesis-generation` | 属立项期（文献搜索、5W1H、Gap Analysis、研究问题卡片）。本课题已进入 P2–P8 执行期，暂不需要。 |
| `rigor-reviewer`（ara-rigor-reviewer） | 要求 ARA 目录结构（PAPER.md、logic/claims.md、trace/exploration_tree.yaml），本项目无此结构；只借鉴其六维尺度中的 D1/D3。 |
| `mathmodel-figure-templates` | 模板固定（SHAP、ROC、泰勒图、环形热图、和弦图等），且使用确定性模拟数据，与本课题波形/扫描曲线类型不符。 |
| `nature-polishing` 的 LaTeX 排版路由 | 需要 `.tex` 工程；本项目交付 docx/pptx。其散文润色部分可用，排版部分不用。 |
| `academic-plotting`（codex） | 架构图走 Gemini 图像生成，非可编辑矢量；数据图能力已被 `nature-figure` + `matplotlib` 覆盖，且其目标为 NeurIPS/ICML 论文。 |
| `3coding-visual` / `2analysis-modeling` / `1start-mathmodel` / `5writing` | 数学建模竞赛的段落式流水线，产物契约（`reports/RESULTS_REPORT.md`、`figures/`）与本研究目录不兼容。 |
| `red-history-fact-checking` | 与本题无关。 |
| `academic-research-suite` / `nature-academic-search` / `nature-reviewer` 完整三审流程 | 面向期刊投稿（三位盲审 + 编辑视角）。答辩阶段只需单审视角，完整流程成本过高，建议只用其 forensic consistency audit 子项。 |

## 5. 补充说明

- `~/.codex/skills/` 下 `scientific-visualization`、`peer-review`、`scientific-writing`、`venue-templates` 四个目录**无 SKILL.md**，仅有 `references/`/`assets/`，故未纳入映射。
- `~/.workbuddy/plugins/cache/` 扫描结果：`cb_teams_marketplace`（document-skills、finance-data）、`codebuddy-plugins-official`（agent-browser、find-skills、playwright-cli、skills-sec-audit）、`experts`（小程序）、`workbuddy-builtin`（设计/云服务/办公等），**未见科研专用 skill**。
- codex 版 `nature-figure` 与 workbuddy 版同名但内容不同（codex 版为单文件，workbuddy 版为 manifest + 契约 + QA 脚本分层），以 workbuddy 版为准。
