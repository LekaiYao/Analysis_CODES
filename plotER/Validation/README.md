# MC validation / sPlot

| 目录/文件 | 用途 |
|---|---|
| [macros/](macros/) | ROOT分析与绘图实现 |
| [workflows/](workflows/) | 参数解析、流程执行、导出检查与provenance工具 |
| aux.h | 原有共用变量范围、分箱和绘图定义 |
| splot_provenance_targets.json | provenance检查对象配置 |
| results/、WEIGHTS/、COMPARE/ | 已有产物，保持原路径与metadata |
| logs/、.cache/ | 历史日志与编译缓存 |

## 常用入口

| 任务 | 入口 |
|---|---|
| 原生1D DATA signal vs MC | [run_DataSIGNAL_VS_MC.sh](workflows/run_DataSIGNAL_VS_MC.sh) |
| 原生2D比较 | [run_DataSIGNAL_VS_MC_2D.sh](workflows/run_DataSIGNAL_VS_MC_2D.sh) |
| 质量相关性 / 稳定性 | [run_mass_correlation.sh](workflows/run_mass_correlation.sh)、[run_mass_shape_stability.sh](workflows/run_mass_shape_stability.sh) |
| PbPb ψ(2S) sPlot流程 | [psi2s_pbpb_splot_validation_workflow.py](workflows/psi2s_pbpb_splot_validation_workflow.py) |
| sWeight导出/核验 | workflows/run_export_*.sh、[validate_sweight_tree.py](workflows/validate_sweight_tree.py) |
| provenance汇总 | [summarize_splot_provenance.py](workflows/summarize_splot_provenance.py) |
| ppRef ψ(2S) prompt-enriched | [ppref_psi2s_prompt_fit.py](../../scripts/workflows/ppref_psi2s_prompt_fit.py)、[ppref_psi2s_prompt_plots.py](../../scripts/workflows/ppref_psi2s_prompt_plots.py) |
| ppRef X MC/sPlot/sideband | [plot_ppref_x_mc_splot_sideband.py](../../scripts/workflows/plot_ppref_x_mc_splot_sideband.py) |
| PbPb23 X MC/sideband | [plot_pb23_x_mc_sideband.py](../../scripts/workflows/plot_pb23_x_mc_sideband.py) |

上述shell可从仓库根用完整相对路径运行，会定位回原Validation工作目录。直接调用宏时仍从Validation目录执行 `root ... macros/Name.C(...)`。历史脚本中已有的样本硬编码保持原样；文件被整理到目录中不表示它自动成为nominal。

67项变量物理字典（本地资料：`docs/reference/mc_validation_variable_dictionary.md`，不随Git发布） · 目录迁移对照（本地资料：`docs/guides/fit_plot_script_layout.md`，不随Git发布）
