#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""Delivered-tag empirical Punzi a=b=5, score 0.50--0.95 in steps of 0.01."""
from pathlib import Path
import json,math,csv,hashlib,shutil,subprocess,argparse,re
from decimal import Decimal,localcontext
import numpy as np,uproot,numexpr,ROOT as R
from fit_input_contract import selection_expression,load_training
ap=argparse.ArgumentParser();ap.add_argument('--tag',required=True);args=ap.parse_args()
repo=Path(__file__).resolve().parents[2];tag=args.tag;root=repo/'fitER/results/ml_fits'/tag;out=root/'punzi_corrected_score';assert not out.exists()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
t,registered=load_training(repo,tag);assert t['weight']=='pThatreweight';base=Path(t['selected_dir'])
specs={**t['sources'],'reference':dict(path=str(base/'REFERENCE_MC_with_score.root'),tree='ntmix_X3872')};arrays={};inputs={};weights={}
for key,spec in specs.items():
 path=Path(spec['path']);inputs[key]=dict(path=str(path),sha256=sha(path),size=path.stat().st_size)
 selection=t['selection'] if key=='data' else t.get('mc_apply_selection',t['selection']);cols=sorted(set(['Bmass','Prediction']+re.findall(r'\bB\w+\b',selection)+([] if key=='data' else ['pThatreweight'])))
 with uproot.open(path) as f:raw=f[spec['tree']].arrays(cols,library='np')
 assert all(np.isfinite(v).all() for v in raw.values());assert numexpr.evaluate(selection_expression(selection),local_dict=raw).all();arrays[key]={k:v.astype(float) for k,v in raw.items()}
 if key!='data':
  assert (raw['pThatreweight']>0).all();values,counts=np.unique(raw['pThatreweight'],return_counts=True);weights[key]=[dict(weight=float(w),entries=int(n)) for w,n in zip(values,counts)]
mc,data,ref=[arrays[k] for k in ['mc','data','reference']];md=(mc['Bmass']>3.62)&(mc['Bmass']<4);den=math.fsum(mc['pThatreweight'][md]);assert den>0
order=np.argsort(ref['Prediction'],kind='stable');cum=np.cumsum(ref['pThatreweight'][order])/ref['pThatreweight'].sum()
for point in registered:
 idx=min(int(np.searchsorted(cum,1-point['target_efficiency'],side='left')),len(order)-1);cut=float(ref['Prediction'][order[idx]]);assert cut==point['score_threshold']
 achieved=math.fsum(ref['pThatreweight'][ref['Prediction']>cut])/math.fsum(ref['pThatreweight']);assert abs(achieved-point['achieved_efficiency'])<1e-12
 if 'data_entries' in point:assert int(np.sum((data['Prediction']>cut)&(data['Bmass']>3.62)&(data['Bmass']<4)))==point['data_entries']
contract=dict(tag=tag,training=t,inputs=inputs,selection=t['selection'],mc_apply_selection=t.get('mc_apply_selection',t['selection']),weight='pThatreweight',weight_values=weights,mc_sumw_D_no_score_cut=den,mass_range=[3.62,4.00],window=[3.84,3.90],signal_window_factor=1,sideband_windows=[[3.75,3.84],[3.90,3.99]],background_width_ratio=1/3,comparison='Prediction > threshold; strict mass boundaries',reference_yield_used=False,reference_mean_used=False)
left=(data['Bmass']>3.75)&(data['Bmass']<3.84);right=(data['Bmass']>3.90)&(data['Bmass']<3.99);refsum=math.fsum(ref['pThatreweight']);rows=[]
for i in range(50,96):
 cut=i/100;passed=data['Prediction']>cut;nl=int(np.sum(left&passed));nr=int(np.sum(right&passed));B=(nl+nr)/3;sw=math.fsum(mc['pThatreweight'][md&(mc['Prediction']>cut)]);fraction=sw/den
 D=25/8+225/13+5*math.sqrt(B)+2.5*math.sqrt(25+20*math.sqrt(B)+4*B);fom=fraction/D
 # Independent sorted-mass count and high-precision formula.
 masses=np.sort(data['Bmass'][passed]);counts=[int(np.searchsorted(masses,h,'left')-np.searchsorted(masses,l,'right')) for l,h in [(3.75,3.84),(3.90,3.99)]];assert counts==[nl,nr]
 with localcontext() as ctx:
  ctx.prec=50;q=Decimal(nl+nr)/3;DD=Decimal(25)/8+Decimal(225)/13+5*q.sqrt()+Decimal('2.5')*(25+20*q.sqrt()+4*q).sqrt();assert abs(float(Decimal(str(fraction))/DD)-fom)<1e-15
 rows.append(dict(score_threshold=cut,xeff_actual_percent=100*math.fsum(ref['pThatreweight'][ref['Prediction']>cut])/refsum,left_entries=nl,right_entries=nr,background=B,mc_sumw_D=sw,mc_fraction_D=fraction,Smin_corrected=D,punzi=fom))
best=max(rows,key=lambda r:r['punzi']);summary={**best,'points':46,'above_99percent_peak':[r['score_threshold'] for r in rows if r['punzi']>=.99*best['punzi']],'above_95percent_peak':[r['score_threshold'] for r in rows if r['punzi']>=.95*best['punzi']],'at_grid_boundary':best in [rows[0],rows[-1]]}
out.mkdir(parents=True);art=out/'artifacts';art.mkdir();dump(out/'scan.json',rows);dump(out/'summary.json',summary)
with (out/'scan.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
reference=repo.parent/'AN/figures/selectionER/ntmix_optimalCUT/bdt_optimization_punzi_ppRef_Bpt_7p5_50p0_inclusive.pdf'
contract.update(a=5,b=5,metric='MC weighted fraction / corrected Smin (Punzi Eq.8)',Smin='a^2/8+9*b^2/13+a*sqrt(B)+b/2*sqrt(b^2+4*a*sqrt(B)+4*B)',score_grid=[i/100 for i in range(50,96)],xeff_grid=None,thresholds_sha256=t['thresholds_sha256'],plot_reference=dict(path=str(reference),sha256=sha(reference)),error_bars=False,new_fits=0,default_changed=False);dump(out/'contract.json',contract)
plot_rows=[r for r in rows if r['score_threshold']<0.95]
ymin=.9*min(r['punzi'] for r in plot_rows);ymax=1.1*max(r['punzi'] for r in plot_rows)
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0)
c=R.TCanvas('c_v27_corrected','',800,600);c.SetLeftMargin(.14);frame=c.DrawFrame(.5,ymin,.95,ymax);frame.SetTitle(' ; Prediction; FOM')
g=R.TGraph(len(plot_rows))
for i,r in enumerate(plot_rows):g.SetPoint(i,r['score_threshold'],r['punzi'])
g.SetName('punzi_corrected_score');g.SetMarkerStyle(20);g.SetMarkerColor(R.kAzure+1);g.SetLineColor(R.kAzure+1);g.SetLineWidth(2);g.Draw('LP SAME')
line=R.TLine(best['score_threshold'],ymin,best['score_threshold'],best['punzi']);line.SetLineStyle(2);line.SetLineWidth(2);line.SetLineColor(R.kAzure+1);line.Draw()
label=R.TLatex();label.SetNDC();label.SetTextFont(42);label.SetTextSize(.035);label.DrawLatex(.16,.86,tag.replace('X_pb23_','PbPb23 ').replace('_rw0_xgb_',' '));label.DrawLatex(.16,.81,'Corrected Punzi, a = b = 5');label.DrawLatex(.16,.76,f"Best threshold = {best['score_threshold']:.2f}");c.SaveAs(str(out/'punzi_score.pdf'))
f=R.TFile(str(art/'curve.root'),'RECREATE');g.Write();f.Close();shutil.copy2(Path(__file__),art/Path(__file__).name)
for key,info in contract['inputs'].items():assert sha(info['path'])==info['sha256']
dump(out/'audit.json',dict(status='PASS',points=46,independent_counts_and_decimal_formula=True,registered_thresholds_exact=len(registered),source_hashes_unchanged=True,no_error_bars=True,root_graph_type='TGraph'))
(out/'README.md').write_text(f'''# {tag}：经验修正式score细扫描

[主图：Punzi vs score](punzi_score.pdf) · [逐点CSV](scan.csv) · [最优点/平台](summary.json) · [核验](audit.json)

固定a=b=5，采用Punzi经验修正式（原文式8），FOM=MC加权通过比例/Smin_corrected。score=0.50–0.95含两端、步长0.01，共46点；严格Prediction>cut。

最大点：score>{best['score_threshold']:.2f}，对应实际加权REFERENCE xeff={best['xeff_actual_percent']:.6f}%，FOM={best['punzi']:.10g}；侧带{best['left_entries']}+{best['right_entries']}条，B={best['background']:.6f}。
99%峰高点集：{summary['above_99percent_peak']}；95%峰高点集：{summary['above_95percent_peak']}。这些集合不是置信区间，最大点仅指本次离散网格。

沿用该tag自身DATA/MC筛选及pThatreweight，不额外增删权重档；具体selection和权重类别见contract.json。MC求和范围3.62<m<4.00，信号窗因子1；SB(3.75,3.84)/(3.90,3.99)，背景为SB计数/3。无参考DATA产额或峰位输入，无拟合/bootstrap，未改变默认。每输入一次读取、46点共用数组；source/selection/weight/窗口变化使结果失效，无需新增大型缓存。

绘图只复用AN参考图样式（亮蓝色实点连线、不展示score=0.95的点，纵轴为可见点最小值×0.9至最大值×1.1、最优点竖虚线），不套用AN的物理配置；横轴限定0.50–0.95，不画error bar，不额外归一化。参考文件及hash见contract.json。
''')
print(json.dumps(summary,indent=2))

assert sha(t['thresholds_file'])==t['thresholds_sha256']
assert sha(base/'batch_apply_summary.json')==t['batch_apply_summary_sha256']
