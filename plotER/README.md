# plotER：绘图入口

> 发布范围：本 README 中的 results、ROOTfiles、archive、backups、metadata 产物及 docs 链接为 CERN 本地资料，不随 Git 发布；源码、配置目录的链接随仓库提供。

[当前 Validation](Validation/README.md) · [2026-09-28 旧 pThat 权重阶段归档](../archive/2026-09-28/README.md)

旧Validation、普通绘图、预筛选图、日志已归档。当前 ppRef 与 PbPb23 B DATA 绘图入口见下表。历史绘图模板、分箱、范围、布局和误差表示保留供后续复用。

后续PbPb23/24 B DATA优先使用已确认pre-cut的六份预处理ROOT，输入与新侧带见[接口说明](../docs/guides/interface_commands.md#pbpb-b-data预处理输入2026-10-14后续绘图默认)。旧raw/fiducial图保留历史定义。

| 目录 | 用途 |
|---|---|
| [results/prompt_x_pb23_pb25/](results/prompt_x_pb23_pb25/README.md) | PbPb23/25 prompt X：含/无pThat5的六变量年份对比，单位面积及23/25比值 |
| [results/Bvar_pb23_sbMC/](results/Bvar_pb23_sbMC/README.md) | Bu/Bs/Bd pb23 当前五项基础cut下，新侧带DATA vs MC变量分布 |
| [results/Bvar_pb23_MC/](results/Bvar_pb23_MC/README.md) | Bu/Bs/Bd PbPb23 MC基础fiducial变量分布，pThatreweight加权 |
| [results/Bmass_PbPb23/](results/Bmass_PbPb23/README.md) | Bu/Bs/Bd PbPb23 DATA：基础fiducial及再加2D>4质量谱，灰色ML侧带边界 |
| [results/Bmass_ppRef/](results/Bmass_ppRef/README.md) | Bu/Bs/Bd 原始 ppRef DATA 无额外 cut 的 Bmass baseline |
| [results/Bvar_ppRef_dataMC_2D4/](results/Bvar_ppRef_dataMC_2D4/README.md) | ppRef B fiducial + 2D>4 DATA/MC变量对比：Bu完成，Bs/Bd待新版MC |
| [results/Variables_ppRef/](results/Variables_ppRef/README.md) | Bu/Bs/Bd 原始 ppRef DATA 变量图与合订本 |
| [macros/](macros/) | DATA/MC、质量谱、non-prompt位移绘图宏 |
| [workflows/](workflows/) | 绘图shell入口 |
| [validation_archive/](validation_archive/README.md) | 历史Validation实现及日期归档导航 |
| correlations/、aux/ | 相关性脚本、共用参数与质量常量 |
| [results/pb18_an_score/](results/pb18_an_score/README.md) | 2018 DATA：AN2019 的 Score>0.06 质量谱 |

普通宏仍从plotER工作目录调用。新样本输入与权重由版本化交付确定，不直接复用旧缓存或硬编码历史输入。归档未更改绘图算法或物理定义。

B 绘图的版本化配置及本地输入准备见 [配置说明](configs/b_inputs/README.md)。
