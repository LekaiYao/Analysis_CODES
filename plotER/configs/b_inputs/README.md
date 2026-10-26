# B 绘图输入配置

这里从已验收的本地结果提取变量、分箱、范围与顺序；内容和来源 SHA 见 `sources.json`。不改变 selection、权重、图形布局或误差定义。原始结果及 provenance 保持原样。

B inventory 位于相邻 ML 项目的 `output/internal/code_validation/b_samples_request_20261013/inventory.json`。它与原结果冻结的 inventory 逐字节一致；SHA-256 为 `cfba30f234d23e13c4620b02b2589a0038d1dbb478f629da634923203f0ded02`。仓库不附带 ML inventory、ROOT 或预处理 manifest；使用者需从相应的版本化 ML 交付获取，并核对源码执行的 size/mtime 和 manifest 检查。不能用新样本静默替代历史输入。

`plot_ppref_b_baseline.py --output <新目录>` 生成 baseline。`plot_ppref_b_mass_2d4.py`、`plot_ppref_b_mass_fid_2d4.py` 是已有结果的增量工具：需要本地 `plotER/results/Bmass_ppRef/metadata/{provenance,variables}.json` 和各 channel 的 `pdf/baseline.pdf`；后者还需要 `pdf/2D4.pdf` 和源码指定的历史 ML `run_metadata.json`。这些输入应从同一份已验收 baseline 交付恢复，不能拼接不同版本。脚本仍使用原暂存目录，执行前应核对路径；本次整理不运行这些工具。

PbPb DATA 新侧带绘图通过 `b_pbpb_preprocessed_input.py` 读取相邻项目 `input/B_pbpb_precut_data.manifest.json`。raw/fiducial 历史脚本仍保留原选择，不因路径修正升级为新流程。
