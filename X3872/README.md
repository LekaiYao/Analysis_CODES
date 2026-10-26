# X3872：去除来源5贡献后的训练输入

> 发布范围：本 README 中的 results、ROOTfiles、archive、backups、metadata 产物及 docs 链接为 CERN 本地资料，不随 Git 发布；源码、配置目录的链接随仓库提供。

2026-10-22：[PbPb23/PbPb25三权重交付](three_weights_pb23_pb25/README.md)已8/8文件、24树PASS；按用户批准从当前五档权重反演L并检验版本内自洽。原权重及旧四份PbPb23 MC保持不变，另有[完整备份](backups/pb23_before_three_weights/backup_manifest.json)，未经用户确认不删除。下文仍为2026-09-29历史交付定义。

**2026-09-29交付完成：8份MC、2份DATA；全部核验通过。** 来源5条目全部保留且权重为0。目录与flatER平行，文件名沿用Henrique的命名。

- 2023：[PbPb23/](PbPb23/)
- 2018：[legacy_run2/REflated/](legacy_run2/REflated/)

每年份包含prompt/nonprompt X及ψ(2S)四份MC。DATA为原文件逐字节复制，SHA一致。

## 权重规则与来源

来源字段为`pThats`，实际产生值为`pthat`，输出权重为`pThatreweight`。RECO和GEN采用同一规则：

- `pThats == 5`或`pthat < pThats`：权重0，不删除条目。
- 其余条目按实际pthat落入[10,15)、[15,30)、[30,50)、[50,∞)，用该区间覆盖的来源10/15/30/50累计贡献计算权重。恰好等于来源下限有效。
- 2023使用Henrique提交4855f0156052a15096e04c07acb466e658940179的[归一化表](https://github.com/Legoinha/Analysis_CODES/blob/4855f0156052a15096e04c07acb466e658940179/flatER/config/mc_normalization.csv)。每来源贡献`L = N_gen / (sigma * filter_eff)`，每区间`weight = 1 / sum(L)`。N_gen使用表中的完整生成数，不用选后候选数。
- 2018按用户确认的假设，从原五档正权重反推累计贡献；每个保留区间`w_new = 1 / (1/w_old - 1/w5)`，保留原总体归一化。此反推不等于找回2018原始截面/效率生产表。

2023只修改weight。2018从原始事件恢复`pthat`(Float_t)与`pThats`(Int_t)，同时修复旧GEN全1权重；其余变量逐值复制。两年份树名、标题、原分支类型、条目数和顺序均保持。

## 核验

| MC | RECO条目 | RECO零权重 | GEN条目 | GEN零权重 |
|---|---:|---:|---:|---:|
| 2023 prompt X | 252,044 | 3,486 | 2,535,361 | 503,859 |
| 2023 nonprompt X | 210,463 | 2,502 | 2,312,537 | 304,363 |
| 2023 prompt ψ(2S) | 202,128 | 1,911 | 2,219,691 | 314,069 |
| 2023 nonprompt ψ(2S) | 168,771 | 1,528 | 2,205,214 | 280,659 |
| 2018 prompt X | 285,732 | 37,902 | 2,524,410 | 515,184 |
| 2018 nonprompt X | 158,490 | 1,631 | 2,191,535 | 289,750 |
| 2018 prompt ψ(2S) | 167,215 | 1,461 | 2,253,734 | 331,925 |
| 2018 nonprompt ψ(2S) | 139,313 | 1,115 | 2,297,386 | 291,622 |

16棵树全部逐列核对通过；独立高精度公式重算所有权重通过；ROOT独立检查零权重和指定fiducial计数通过。输入/输出SHA、来源表及脚本版本保存在[构建记录](metadata/mc_manifest.json)，复核见[独立验收](metadata/independent_acceptance.json)。零权重列包含来源5及其它来源低于自身阈值的条目。

2018四种MC的RECO全部唯一匹配原始候选，旧权重逐条一致；GEN全部原列和顺序与既有flatten选择一致，prompt ψ另以PVxyz定位独立复核。已有映射缓存避免重复读取大型原始文件；输入版本或flatten选择变化即失效，无需另建缓存。

## ML使用

**除修改输入路径，需在MC训练、APPLY/参考和Punzi统一使用`pThatreweight > 0`；DATA不加此条件。** 当前训练接受零权重，但Punzi拒绝零权重。

2018当前三变量`BtrkLeadPt/BtrkMaxdR/BtrkSubPt`均可用；还需登记2018年份/路径。2018缺`nChargedTracks_LOOSE/TIGHT`，旧`nChargedTracks`为int32常数100，旧`Score/mvapref`照原样保留，不能将旧Score作为新模型Prediction。

生成命令（已有输出受防覆盖保护；复现需另设空输出目录）：

```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-opt/setup.sh
python3 X3872/build_reweighted_samples.py --source5-policy zero --output X3872 --maps X3872/metadata/run2_source_maps --writer X3872/ReweightFlatTree.C
python3 X3872/audit_reweighted_samples.py --root X3872
```

未修改上游ROOT或ML配置，未训练。详细权重历史见[权重说明](../docs/reference/pthat_reweighting.md)。

## 源码复现前提

`build_reweighted_samples.py` 是 2026-09-29 的历史 2018/2023 流程；默认 writer 为仓库的 `ReweightFlatTree.C`，默认映射为本地 `metadata/run2_source_maps/`。需要恢复每种 promptness/particle 的 `reco.npz`、`gen.npz`、`acceptance.json`，以及本地 `docs/reference/upstream_review_evidence/2026-09-29/pthat_layers/pthat_layers_audit.json`。这些来自同一次已验收交付；映射通过源文件 size/mtime 检查，不可跨版本复用。历史 CSV 已版本化，脚本从原提交 `4855f0156052a15096e04c07acb466e658940179` 读取，因此 clone 还需包含该 Git 对象。缺少映射时，`metadata/prepare_2018_weight_mapping.py` 是历史生成工具，输出仍为 `/tmp/pthat_rebuild_maps`；它需要完整原始样本，不能把生成映射当作轻量检查。

`build_three_weights.py` 是 2026-10-22 的两阶段记录式生产工具。输入为本地 `three_weights_pb23_pb25/metadata/input_inventory.json`、同目录历史 CSV 和 inventory 指定的 backup manifest；从完整已验收交付恢复这些输入及其 SHA 对应的 ROOT。它们不随 Git 上传。第一阶段只处理已明确归一化的样本，留下 `PARTIAL_BLOCKED_NORMALIZATION`；`--resume-from-henrique-bands` 只用于这一状态及已有三个产物，在当时确认的反演假设下完成其余五份。已完成的 PASS 目录不能重新运行 resume。

这不是复制 inventory 即可对任意新目录生产的 CLI：inventory 含绝对输出路径和旧输入指纹。新生产需先准备并审核新 inventory、独立输出及备份，不能把历史输出路径直接复制执行。当前 PbPb23/25 对比入口为 `scripts/workflows/plot_prompt_x_pb23_pb25.py`；旧 PbPb18/23 脚本保留历史定义。本轮仅发布源码，未生成权重、未执行拟合、未改动科研结果。
