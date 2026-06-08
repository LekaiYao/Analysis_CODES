# Bootstrap：当前优化、测试与拟合工作流

## 范围
本文档描述当前 `selectionER/` 与 `fitER/` 中实际实现、并已在本仓库内落地的工作流。

主要组件：
- `selectionER/optimalCUT_punzi.C`
- `selectionER/optimalCUT_punzi_test.C`
- `selectionER/punzi_test_matrix.conf`
- `selectionER/optimalCUT_fom.C`
- `selectionER/optimalCUT.conf`
- `selectionER/run_punzi_test_condor.sh`
- `selectionER/submit_punzi_test_condor.sh`
- `fitER/bmeson_fit_from_conf.sh`
- `fitER/roofitB.C`

## 1) 共享配置
主配置文件：
- `selectionER/optimalCUT.conf`

结构：
- 多 profile section；每个 tag 对应一个样本或一套配置。

典型键值：
```ini
[Bs_pp24_v1_fid1_17v1_xgb_v1]
system=pp
dataPath=...
mcPath=...
dataTreeName=ntphi
mcTreeName=ntphi
scoreVar=xgb_score
preCut=(...)
sidebandLow=(...)
sidebandHigh=(...)
fsRegion=(...)
refScoreCut=0.6
signalWidth=0.16
sidebandWidth=0.2
punziA=2.0
punziB=5.0
fileNamePattern=punzi_....
mass_range=(Bmass > 5.1 && Bmass < 5.7)
bin_width=0.005
optimalCUT_punzi=0.910000
optimalCUT_fom=0.720000
```

当前行为：
- `optimalCUT_punzi.C` 会更新或插入 `optimalCUT_punzi`
- `optimalCUT_fom.C` 会更新或插入 `optimalCUT_fom`
- `optimalCUT_punzi_test.C` 不回写 `optimalCUT.conf`，只输出测试结果

## 2) Punzi 主优化
脚本：
- `selectionER/optimalCUT_punzi.C`

运行方式：
```bash
root -l -b -q 'selectionER/optimalCUT_punzi.C("selectionER/optimalCUT.conf","<profile>")'
```

所需配置键：
- `dataPath`, `mcPath`
- `dataTreeName`, `mcTreeName`
- `scoreVar`
- `preCut`
- `sidebandLow`, `sidebandHigh`
- `signalWidth`, `sidebandWidth`
- `punziA`, `punziB`
- `outputDir`, `fileNamePattern`
- `system`

当前实现定义：
- `sTotal = N_MC(preCut)`
- 对每个阈值 `thr`，从 `0.00` 到 `0.99`，步长 `0.01`：
- `s = N_MC(preCut && scoreVar > thr)`
- `bLow = N_DATA(sidebandLow && preCut && scoreVar > thr)`
- `bHigh = N_DATA(sidebandHigh && preCut && scoreVar > thr)`
- `sigEff = s / sTotal`
- `sbToSigScale = signalWidth / sidebandWidth`
- `bkg = (bLow + bHigh) * sbToSigScale`
- `S_min = punziSmin(bkg, punziA, punziB)`
- `Punzi_FOM = S_min / sigEff`

当前 `punziSmin` 公式：
```text
S_min = b^2 / 2 + a * sqrt(B) + b / 2 * sqrt(b^2 + 4 a sqrt(B) + 4 B)
```
其中：
- `a = punziA`
- `b = punziB`
- `B = bkg`

选择规则：
- 最小化 `Punzi_FOM`
- 若 `sigEff <= 0`，跳过该点
- 不扫描 `thr = 1.0`

输出：
- 图保存到 `outputDir/fileNamePattern`
- 最优阈值会回写到同一 profile 的 `optimalCUT_punzi=<value>`

## 3) Punzi Test Workflow
脚本：
- `selectionER/optimalCUT_punzi_test.C`

用途：
- 在不改动主 `optimalCUT_punzi.C` 的前提下，对多套 `S_min` 公式和多组 `(a,b)` 参数做批量扫描测试
- 不回写 `optimalCUT.conf`
- 输出图和日志到测试目录

### 3.1 本地串行运行
从 `selectionER/` 目录运行：
```bash
cd selectionER
root -l -b -q 'optimalCUT_punzi_test.C("optimalCUT.conf","<profile>")'
```

示例：
```bash
cd selectionER
root -l -b -q 'optimalCUT_punzi_test.C("optimalCUT.conf","Bd_pp24_v1_fid1_9v1_xgb_v1")'
```

当前实现定义：
- 计算方式仍沿用主 Punzi workflow：
- `sTotal = N_MC(preCut)`
- `bkg = (bLow + bHigh) * signalWidth / sidebandWidth`
- `Punzi_FOM = S_min / sigEff`
- 扫描阈值 `thr = 0.000, 0.002, ..., 0.998`
- 当前图只画曲线和最优 cut 竖线，不再描点

### 3.2 外置参数矩阵
参数矩阵文件：
- `selectionER/punzi_test_matrix.conf`

格式：
```text
# formula,file_tag,a,b,enabled
simplified,simplified,1.64,3.0,1
gaussian,gaussian,1.64,5.0,1
improved,improved,2.0,5.0,1
```

字段含义：
- 第 1 列：公式名，只能是 `simplified`、`gaussian`、`improved`
- 第 2 列：输出文件标签 `file_tag`
- 第 3、4 列：`a`、`b`
- 第 5 列：`enabled`
  - `1` 表示参与运行
  - `0` 表示跳过

当前支持的 `S_min` 公式：
- `simplified`
```text
S_min = a/2 + sqrt(B)
```
- `gaussian`
```text
S_min = b^2/2 + a sqrt(B) + b/2 * sqrt(b^2 + 4 a sqrt(B) + 4 B)
```
- `improved`
```text
S_min = a^2/8 + 9 b^2/13 + a sqrt(B) + b/2 * sqrt(b^2 + 4 a sqrt(B) + 4 B)
```

### 3.3 单配置筛选运行
测试宏现在支持额外参数：
```bash
cd selectionER
root -l -b -q 'optimalCUT_punzi_test.C("optimalCUT.conf","<profile>","punzi_test_matrix.conf","<selection_key>")'
```

其中 `selection_key` 规则是：
```text
<file_tag>_a<aa>_b<bb>
```
例如：
- `simplified_a1p64_b3p00`
- `gaussian_a1p64_b5p00`
- `improved_a2p00_b5p00`

用途：
- 本地只跑一组配置
- Condor 中每个 job 只跑一组配置，实现并行

### 3.4 输出目录结构
若不传 `selection_key`，输出目录为：
```text
selectionER/opt_tests/punzi_scan/<profile>/
```
其中包含：
- 多张 PDF
- `scan.log`
- `summary.log`

若传入 `selection_key`，输出目录为：
```text
selectionER/opt_tests/punzi_scan/<profile>/<selection_key>/
```
其中包含：
- 单张 PDF
- `scan.log`
- `summary.log`

## 4) Punzi Test Condor Workflow
### 4.1 设计目标
当前仓库代码与输出路径在 EOS 下；但 lxplus 上进行 Condor 提交时，实际提交目录通常应放在 AFS 下。

当前实现采用“两地分离”的方式：
- `repo_root` 指向 EOS 仓库
- `submit dir` 建在 AFS
- Condor worker 实际回到 EOS 仓库里的 `selectionER/` 执行 ROOT 宏

### 4.2 Worker Wrapper
脚本：
- `selectionER/run_punzi_test_condor.sh`

用途：
- 每个 Condor job 只跑一个 `selection_key`
- 支持以下参数：
```bash
run_punzi_test_condor.sh <repo_root> <profile> <selection_key> [conf_path] [matrix_path]
```

若 worker 节点上 `root` 不在 PATH 中，可通过环境变量：
- `ROOT_SETUP_SCRIPT`
来提供 ROOT 环境初始化脚本。

### 4.3 Submit Helper
脚本：
- `selectionER/submit_punzi_test_condor.sh`

用途：
- 从 `punzi_test_matrix.conf` 中展开所有 `enabled=1` 的配置
- 为每个配置生成一个 `selection_key`
- 在 AFS 下生成 submit 目录、jobs 清单、wrapper 副本和 Condor submit 文件

运行方式：
```bash
cd selectionER
bash submit_punzi_test_condor.sh <profile> <afs_submit_root> [repo_root] [conf_path] [matrix_path]
```

示例：
```bash
cd selectionER
bash submit_punzi_test_condor.sh Bd_pp24_v1_fid1_9v1_xgb_v1 /afs/cern.ch/user/l/leyao/private/punzi_condor
```

当前行为：
- 要求 `afs_submit_root` 必须是 `/afs/...` 路径
- 在 AFS 下创建一次 submit run 目录
- 生成：
  - `jobs.txt`
  - `submit_punzi_test.sub`
  - `run_punzi_test_condor.sh`
  - `condor_logs/`
- submit 文件采用：
  - `queue selection_key from jobs.txt`
- 每个 `selection_key` 形成一个独立 Condor job

脚本最后会打印推荐提交步骤：
```bash
cd <generated_run_dir>
module load lxbatch/eossubmit
condor_submit submit_punzi_test.sub
```

## 5) FOM 优化
脚本：
- `selectionER/optimalCUT_fom.C`

运行方式：
```bash
root -l -b -q 'selectionER/optimalCUT_fom.C("selectionER/optimalCUT.conf","<profile>")'
```

所需配置键：
- `dataPath`, `mcPath`
- `dataTreeName`, `mcTreeName`
- `scoreVar`
- `preCut`
- `sidebandLow`, `sidebandHigh`
- `fsRegion`
- `refScoreCut`
- `fileNamePattern`
- `system`

当前定义：
- `sideband = sidebandLow || sidebandHigh`
- `fsWidth` 从 `fsRegion` 解析
- `sidebandWidth` 从 `sidebandLow + sidebandHigh` 解析
- `bkgScale = fsWidth / sidebandWidth`

参考归一化：
- `N_MC(ref) = N_MC(preCut && fsRegion && scoreVar > refScoreCut)`
- `N_DATA(ref) = N_DATA(preCut && fsRegion && scoreVar > refScoreCut)`
- `B_refSB = N_DATA(sideband && preCut && scoreVar > refScoreCut)`
- `Fs = N_MC(ref) / (N_DATA(ref) - bkgScale * B_refSB)`

逐阈值扫描：
- 对每个 `thr`，从 `0.00` 到 `1.00`，步长 `0.01`：
- `s = N_MC(preCut && scoreVar > thr)`
- `b = N_DATA(sideband && preCut && scoreVar > thr)`
- `den = Fs * s + bkgScale * b`
- `FOM = (Fs * s) / sqrt(den)`

选择规则：
- 最大化 `FOM`

输出：
- 图保存到 `./opt_plots/`，文件名由 `fileNamePattern` 中的 `punzi` 替换为 `fom`
- 最优阈值会回写到同一 profile 的 `optimalCUT_fom=<value>`

## 6) 从配置驱动 B 介子拟合
脚本：
- `fitER/bmeson_fit_from_conf.sh`

用途：
- 直接从 `selectionER/optimalCUT.conf` 的单个 profile 运行 `Bu`、`Bs` 或 `Bd` 的 FULL 拟合

当前正确运行方式：
```bash
cd fitER
bash bmeson_fit_from_conf.sh <profile>
```

示例：
```bash
cd fitER
bash bmeson_fit_from_conf.sh Bd_pp24_v1_fid1_9v1_xgb_v1
```

当前实现特点：
- `optimalCUT.conf` 路径硬编码为 `../selectionER/optimalCUT.conf`
- 运行前必须在 `fitER/` 目录下
- 不再接收 `conf_path` 参数

所需配置键：
- `dataPath`, `mcPath`
- `dataTreeName`, `mcTreeName`
- `system`
- `scoreVar`
- `preCut`
- `optimalCUT_punzi`
- `mass_range`，格式需为 `(Bmass > a && Bmass < b)`
- `bin_width`

可选键：
- `channel=Bu/Bs/Bd`

若缺少 `channel`，则由 `mcTreeName` 自动推断：
- `ntKp -> Bu`
- `ntphi -> Bs`
- `ntKstar -> Bd`

当前拟合 cut：
```text
(preCut) && (scoreVar > optimalCUT_punzi)
```

重要说明：
- 当前拟合 cut 不再自动注入 sideband cut

质量窗口与 binning：
- shell 脚本会解析 `mass_range`
- 脚本会导出：
  - `ROOFIT_MASS_MIN`
  - `ROOFIT_MASS_MAX`
  - `ROOFIT_BIN_WIDTH`
- `roofitB.C` 会读取这些环境变量并覆盖：
  - `Bmass` 的拟合/作图范围
  - `m_rangeB`
  - `fitMassBins = round((max - min) / bin_width)`

底层 ROOT 调用：
```bash
root -b -q 'roofitB.C++("<TREE>", 1, "<DATA>", "<MC>", "Bpt", "<CUTS>", "<SYSTEM>")'
```

## 7) 当前拟合后端
文件：
- `fitER/roofitB.C`
- `fitER/roofitB.h`
- `fitER/aux/uti.h`

当前状态：
- B 介子拟合路径已经通过配置驱动 wrapper 接通
- 最近可工作的 ROOT 环境为：
  - ROOT `6.32.02`
- `fitER/aux/uti.h` 中包含 `#include "RooChi2Var.h"`

## 8) 已知注意事项
- `optimalCUT_fom.C` 的输出目录当前是硬编码 `./opt_plots`，不是读 `conf` 的 `outputDir`
- `optimalCUT.conf` 里当前存在重复 tag：
  - `[Bs_pb24_v2_fid1_17v1_xgb_v1]`
- 一些 profile 的 `sidebandLow` 是零宽区间，例如：
  - `Bs_pp24_v1_fid1_17v1_xgb_v1`
  - `Bd_pp24_v1_fid1_9v1_xgb_v1`
- `optimalCUT.conf` 中的 ROOT 文件相对路径依赖运行目录；当前 Punzi/FOM/test workflow 应在 `selectionER/` 下运行，B fit workflow 应在 `fitER/` 下运行
