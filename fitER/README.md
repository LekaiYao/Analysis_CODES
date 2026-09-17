# fitER：拟合入口

默认2023 X拟合从仓库根运行：

```bash
python3 scripts/workflows/run_pb23_default_fit.py --tag TAG --output OUTPUT
```

参数与环境见默认workflow（本地资料：`docs/pb23_default_fit_workflow.md`，不随Git发布）。只运行用户指定tag；六个xeff点的模型与扫描定义保持不变。

| 目录/文件 | 用途 |
|---|---|
| [workflows/](workflows/) | Python流程与shell入口；按下表选择 |
| [models/](models/) | 单年份及跨年份X/ψ(2S)拟合模型 |
| [cache/](cache/) | Prepare*输入筛选缓存宏；不是缓存产物目录 |
| [diagnostics/](diagnostics/) | MC峰形诊断、重画、结果导出及历史注入研究工具 |
| roofitB.C、roofitB.h、aux/ | 保留的通用原生拟合核心与共用头文件 |
| configs/ | 已有版本化物理配置 |
| results/（本地资料：`fitER/results/README.md`，不随Git发布）、ROOTfiles/ | 已有输出，路径不变 |
| logs/、.cache/ | 历史运行日志及编译/字节码缓存，非分析入口 |

## 按任务选择

| 任务 | 入口 |
|---|---|
| 当前2023 X默认六点 | [run_pb23_default_fit.py](../scripts/workflows/run_pb23_default_fit.py) |
| X MC定形单年份/兼容扫描 | [x_fit_scan_workflow.py](workflows/x_fit_scan_workflow.py) |
| X两年份MC定形 | [x_mc_shape_simultaneous_year_fit_workflow.py](workflows/x_mc_shape_simultaneous_year_fit_workflow.py) |
| X单年、混合与同时拟合比较 | [run_x_two_year_fit_strategy_comparison.py](workflows/run_x_two_year_fit_strategy_comparison.py) |
| ψ(2S)单年份MC定形 | [psi2s_fit_scan_workflow.py](workflows/psi2s_fit_scan_workflow.py) |
| ψ(2S)两年份MC定形 | [psi2s_simultaneous_year_fit_workflow.py](workflows/psi2s_simultaneous_year_fit_workflow.py) |
| 原生ppRef X / ψ(2S) | [X3872doRoofit.sh](workflows/X3872doRoofit.sh)、[Psi2SdoRoofit.sh](workflows/Psi2SdoRoofit.sh) |
| 原生non-prompt | workflows/*doRoofit_nonPrompt.sh |
| 历史H004/H010/H011/H012与cut扫描 | workflows/run_*；保留其原样本、cut及科研授权要求 |

移动后的shell仍将相对输入/输出定位于原fitER目录。Python流程的命令行参数保持不变；直接运行模型宏时，仍从fitER工作目录调用其新路径，例如 `root -l -b -q 'models/PbPbXEfficiencyFit.C(...)'`，实际参数必须遵循对应workflow。

提交入口仍在scripts/，已经同步新路径。目录整理不授权自动执行历史toys或batch作业。旧命令只需按完整迁移对照（本地资料：`docs/guides/fit_plot_script_layout.md`，不随Git发布）替换脚本路径；不保留同名软链接或多份源码。
