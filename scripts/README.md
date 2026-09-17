# Analysis_CODES 脚本入口

现有脚本暂不批量移动，以保护既有 request、Condor 与复现命令。新脚本必须遵守
`docs/guides/artifact_governance.md`（本地资料：`docs/guides/artifact_governance.md`，不随Git发布）：workflow 进入
`scripts/workflows/`，validation 进入 `scripts/validation/`，一次性调查进入
`scripts/dev/<study_id>/`。

## 当前 X 接口

| 脚本 | 角色 |
|---|---|
| `submit_x_mc_shape_simultaneous_year_fit_manifest.py` | 当前 schema-v2 两年份 X simultaneous MC-shape consumer |
| `submit_x_two_year_fit_strategy_comparison.py` | simultaneous 主结果的 merged/independent cross-check 提交与 preflight |
| `submit_x_nominal_fit_manifest.py` | PbPb24 schema-v2 consumer；DATA-only 只能显式 compatibility |
| `import_ml_manifest.py` | 通用单-tag importer；非 paired-year nominal 主桥梁 |
| `run_x_analysis_from_manifest.py` | 隔离 prototype/compatibility 入口 |

实际两年份 comparison runner 位于 `fitER/workflows/run_x_two_year_fit_strategy_comparison.py`。

## 10v1 混合 Punzi 研究

- `workflows/scan_10v1_merged_punzi_fine.py`：完整reference MC重建共同效率阈值，复验旧点后执行10–40%每1%细扫描与峰点混合拟合；最终报告需结合多初值审计，参见正式研究说明。
- `workflows/scan_10v1_merged_punzi.py`：冻结 cache/threshold 的 `preflight → scan → fit`；用共同效率的混合 sideband Punzi 预选候选，再运行固定信号比例的混合质量谱拟合。科学参数显式传入，拒绝覆盖已有结果；定义与结果见 `docs/archive/previous_workflows/10v1_merged_punzi_optimization.md`（本地资料：`docs/archive/previous_workflows/10v1_merged_punzi_optimization.md`，不随Git发布）。

## ppRef snapshot v1 validation

- `workflows/ppref_x_snapshot_v1_workflow.py`：冻结 ppRef X DATA/MC 的 no-ML mass fit、signed sPlot、
  67 common-scalar empirical-CDF validation 与 protected-path postflight；按
  `preflight → fit → splot → cdf → postflight` 分阶段运行。`cdf` 只复用通过 event alignment
  和 weight closure 的版本化 event-level ROOT；输入、selection、workspace/model 或 schema 变化时
  必须新建 revision，不得复用。

## 治理与 validation 工具

- `validation/build_artifact_migration_inventory.py`：由冻结 spec 生成 metadata-only
  inventory、CSV path map 和 validation；no-overwrite，不移动或删除源 artifact。
- `validation/execute_artifact_migration.py`：只消费已通过的冻结 inventory；迁移前后核对
  metadata，引用中的结果目录保留相对 compatibility symlink，失败时回滚且不删除 artifact；
  post-link validation 还会检查归档结果内的 nested symlink 是否断链。

## Psi2S reference/control

- `workflows/ppref_psi2s_pthat_validation.py`：ppRef新版RECO-only、pThat加权MC定形、DATA sPlot和67变量验证；输出目录存在即拒绝覆盖。


| 脚本 | 角色 |
|---|---|
| `submit_psi2s_simultaneous_year_fit_manifest.py` | 两年份 Psi2S simultaneous consumer |
| `submit_psi2s_pbpb_splot_validation.py` | Psi2S sPlot validation |
| `run_psi2s_closure_pilot.py` | 12/22-variable closure；旧样本 reference |
| `run_psi2s_year_mc_comparison.py` | 旧 Psi2S 年份 MC shape 对照 |
| `run_x_baseline_year_mc_comparison.py` | 旧 X baseline 年份 MC shape 对照；不用于新 revision |

## Historical / compatibility / one-off

以下入口保留复现价值，但不得自动作为新 nominal 工作流：

- `submit_x_simultaneous_year_fit_manifest.py`：schema-v1 DATA-only compatibility；
- `submit_x_ablation_nominal_workflow.py`：已收口的消融扫描；
- `merge_x_punzi_results.py`：旧 Punzi 结果合并；
- `run_psi2s_combined_closure.py`：早期 combined closure；
- `submit_psi2s_nominal_fit_manifest.py`：旧 PbPb24 single-year fit；
- `submit_psi2s_data_gaussian_manifest.py`：DATA-Gaussian candidate study；
- `submit_psi2s_high_efficiency_exploratory.py`：高效率 exploratory；
- `submit_psi2s_xeff45_splot_closure.py`：一次性 xeff45 gate。

脚本名本身不代表当前状态；运行前必须从 `docs/README.md` 和对应 manifest/validation 确认角色。

- `workflows/ppref_psi2s_prompt_fit.py`：原 ψ(2S) 预筛选追加 SV–PV significance < 2，重做 MC/Data fit 与 sPlot。
- `workflows/ppref_psi2s_prompt_plots.py`：复用 ROOT 格式输出 67 变量 PDF、对比原 CDF；消费上述筛选缓存。

- `workflows/run_pb23_four_training_scan.py`：四套2023训练六点单年份MC定形拟合，原宏不变；缓存/编译复用及技术中断续跑。

- `workflows/run_pb23_default_fit.py`：当前默认2023 X单年份MC定形六点扫描（xeff=15/20/25/30/35/40%）；--tag指定训练，冻结阈值、按tag原权重，最大sqrt(q0)纳入比较。

- `workflows/plot_ppref_x_mc_splot_sideband.py`：ppRef X普通预筛选，复用sPlot/model，67变量MC/信号sPlot/原始4–8σ侧带比较。

- `workflows/plot_pb23_x_mc_sideband.py`：PbPb23 prompt X MC与指定3.75–3.85/3.9–4.0侧带背景67变量比较，MC pThat加权。

## 运行范围

从仓库根启动脚本，使用项目验证的 ROOT 6.32.02 / LCG_106 环境；输入 ROOT、上游 ML manifests 和本地科研说明不随 Git 发布。

- `run_pb23_two_peak_fit.py` 与 `run_pb23_conditional_fit.py`：仅对明确指定的tag/选点运行；条件检验的 sqrt(q0) 尚未校准。
- `scan_pb23_run2_fom.py`：固定用于 v27 的已批准研究，依赖 `scan_10v1_merged_punzi.py` 中的侧带工具和现存50%参考结果，拒绝覆盖已有输出。`plot_pb23_run2_fom.py`只读取已冻结扫描。
- `run_pb23_four_training_scan.py`、`run_two_peak_fixed_width_check.py`及10v1脚本保留历史复现用途；不要据文件存在自动启动分析。
