# plotER：绘图入口

| 目录 | 用途 |
|---|---|
| [Validation/](Validation/README.md) | sPlot、MC validation、CDF与closure |
| [macros/](macros/) | 普通DATA/MC、质量谱、non-prompt位移绘图宏 |
| [workflows/](workflows/) | 普通绘图shell入口 |
| correlations/ | 相关性Python脚本 |
| aux/ | 共用绘图参数与质量常量 |
| results/、local_test_outputs/ | 已有输出，未搬动 |
| logs/ | 历史运行日志 |

普通宏从plotER目录调用，例如 `root -l -b -q 'macros/plot_dataMC.C("ntmix","ppRef")'`。保持原有工作目录，避免相对输出路径变化。

[PbPb ψ(2S) selection质量扫描入口](workflows/run_pbpb24_psi2s_selection_mass_scan.sh)。当前ppRef/PbPb MC对比的Python流程仍在 [scripts/workflows/](../scripts/workflows/)，具体入口见Validation说明。绘图格式、物理cut及产物格式遵守原workflow。
