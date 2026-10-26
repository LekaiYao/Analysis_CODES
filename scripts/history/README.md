# 历史诊断源码

这里保留原先仅存在于本地结果 metadata 的独有源码，`source_index.json` 记录原位置与 SHA-256，源码逐字节保留。原结果目录及冻结脚本未移动、未删除。

这些文件是历史研究记录，不是当前 production 入口。可能直接执行、修改既有产物，依赖当时的绝对路径、暂存文件、旧结果和本地参考资料；部分脚本依赖当时的工作流函数版本，不能直接与当前版本混用。重现前必须明确输入版本、执行位置、依赖版本及独立输出，并取得该研究的执行授权；不要通过 import 批量检查它们。

三份 `initial_plot_ppref_b_baseline.py` 是已被 `scripts/workflows/plot_ppref_b_baseline.py` 替代的初始快照，仍随本地结果保存，不另发布重复源码。外部 Bfinder/Bntuple 源码摘录仍随本地 provenance 保存。
