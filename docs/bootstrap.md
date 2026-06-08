# Bootstrap：当前优化与拟合工作流

## 范围
本文档描述当前在 `selectionER/` 与 `fitER/` 中实际实现的工作流。

主要组件：
- `selectionER/optimalCUT_punzi.C`
- `selectionER/optimalCUT_fom.C`
- `selectionER/optimalCUT.conf`
- `fitER/bmeson_fit_from_conf.sh`
- `fitER/roofitB.C`

## 1) 共享配置
配置文件：
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

## 2) Punzi 优化
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

当前定义：
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

选择规则：
- 最小化 `Punzi_FOM`
- 若 `sigEff <= 0`，跳过该点
- 不扫描 `thr = 1.0`

输出：
- 图保存到 `outputDir/fileNamePattern`
- 最优阈值会回写到同一 profile 的 `optimalCUT_punzi=<value>`

## 3) FOM 优化
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

## 4) 从配置驱动 B 介子拟合
脚本：
- `fitER/bmeson_fit_from_conf.sh`

用途：
- 直接从 `selectionER/optimalCUT.conf` 的单个 profile 运行 `Bu`、`Bs` 或 `Bd` 的 FULL 拟合

运行方式：
```bash
bash fitER/bmeson_fit_from_conf.sh <profile> [conf_path]
```

重要运行前提：
- 该脚本内部调用的是 `root -b -q "roofitB.C++(...)"`，因此当前工作目录必须能直接找到 `roofitB.C`
- 若从仓库根目录运行，应先进入 `fitER/`，或后续把脚本改为使用显式路径

当前可工作的调用示例：
```bash
cd fitER
bash bmeson_fit_from_conf.sh Bd_pp24_v1_fid1_9v1_xgb_v1 ../selectionER/optimalCUT.conf
```

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

## 5) 当前拟合后端
文件：
- `fitER/roofitB.C`
- `fitER/roofitB.h`
- `fitER/aux/uti.h`

当前状态：
- B 介子拟合路径已经通过配置驱动 wrapper 接通
- 最近可工作的 ROOT 环境为：
- ROOT `6.32.02`
- `fitER/aux/uti.h` 中包含 `#include "RooChi2Var.h"`

最近测试过的 B 例子：
- Profile：`Bs_pp24_v1_fid1_17v1_xgb_v1`
- Fit cut：
- `((abs(By) < 2.4) && (Bpt > 7.5) && (Bnorm_svpvDistance_2D>2) && (BtrkPtimb<0.3) && (Bchi2Prob>0.02)) && (xgb_score > 0.910000)`
- 质量范围：
- `5.1` 到 `5.7`
- bin 宽度：
- `0.005`

生成输出示例：
- `results/pp/ntphi/Bpt/data_pp_Bpt_5_60_ntphi.pdf`
- `results/pp/ntphi/Bpt/mc_pp_Bpt_5_60_ntphi.pdf`

## 6) Recommended Session Entry Points
Optimization:
```bash
root -l -b -q 'selectionER/optimalCUT_punzi.C("selectionER/optimalCUT.conf","<profile>")'
root -l -b -q 'selectionER/optimalCUT_fom.C("selectionER/optimalCUT.conf","<profile>")'
```

Fit:
```bash
cd fitER
bash bmeson_fit_from_conf.sh <profile> ../selectionER/optimalCUT.conf
bash fitER/bmeson_fit_from_conf.sh <profile>
```
