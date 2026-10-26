# 历史 Validation 实现

> 发布范围：本 README 中的 results、ROOTfiles、archive、backups、metadata 产物及 docs 链接为 CERN 本地资料，不随 Git 发布；源码、配置目录的链接随仓库提供。

旧结果、权重、日志、ROOT和provenance已于2026-09-28归档；这里保留原宏、workflow、aux.h与配置供查阅和模板复用。

- [日期归档总入口](../../archive/2026-09-28/README.md)
- [旧研究图册索引](../../archive/2026-09-28/plotER/validation_archive/archive/previous_validation/README.md)
- [旧ROOT与支持文件](../../archive/2026-09-28/plotER/validation_archive/results/)
- [原路径映射](../../archive/2026-09-28/metadata/path_map.json)
- [保留的宏](macros/) · [保留的workflow](workflows/)

历史实现及冻结脚本中的路径仍保留原执行provenance；历史重现需显式按映射定位输入并指定独立输出。它们不自动成为新样本的默认流程。现有分箱、范围、布局、图例和误差表示应优先复用。
