#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""Frozen count-table fine scan, Corrected Punzi a=5, independent b values."""
from pathlib import Path
import json,math,hashlib,shutil,csv,subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
repo=Path(__file__).resolve().parents[2];tag='X_pb23_v27_fid18_9v9_rw0_xgb_v1';root=repo/'fitER/results/ml_fits'/tag;source=root/'punzi_counting';out=root/'punzi_corrected_b';assert not out.exists()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
inputs={n:sha(source/n) for n in ['xeff_scan.json','contract.json']};rows=json.loads((source/'xeff_scan.json').read_text());assert [r['xeff'] for r in rows]==list(range(15,61));assert all(r['status']=='PASS' for r in rows)
a=5.;bs=[1.645,2.,5.];results=[];summary={}
def denom(B,b):return a*a/8+9*b*b/13+a*math.sqrt(B)+b/2*math.sqrt(b*b+4*a*math.sqrt(B)+4*B)
for b in bs:
 group=[]
 for r in rows:
  B=r['background'];N=r['sideband_entries'];assert abs(B-N/3)<1e-10
  f=r['mc_sumw_D']/json.loads((source/'contract.json').read_text())['mc_sumw_D_no_score_cut'];assert abs(f-r['mc_fraction_D'])<1e-12
  D=denom(B,b);G=D-a*a/8-5*b*b/26;assert abs(G-a*math.sqrt(B)-b*math.sqrt(B+G))<1e-10
  derivative=a/(2*math.sqrt(B))+b*(a/math.sqrt(B)+2)/(2*math.sqrt(b*b+4*a*math.sqrt(B)+4*B));error=f*derivative/D**2*r['background_error']
  step=B*1e-5;numeric=abs((f/denom(B+step,b)-f/denom(B-step,b))/(2*step))*r['background_error'];assert abs(error-numeric)<1e-11
  value=f/D
  assert D>G
  point=dict(a=a,b=b,xeff=r['xeff'],xeff_actual_percent=r['xeff_actual_percent'],score_threshold=r['score_threshold'],background=B,mc_fraction_D=f,Smin=D,punzi=value,background_error=error);group.append(point)
 best=max(group,key=lambda r:r['punzi']);summary[str(b)]={**best,'above_99percent_peak':[r['xeff'] for r in group if r['punzi']>=.99*best['punzi']],'above_95percent_peak':[r['xeff'] for r in group if r['punzi']>=.95*best['punzi']]};results.extend(group)
out.mkdir();art=out/'artifacts';art.mkdir();dump(out/'scan.json',results);dump(out/'summary.json',summary)
with (out/'scan.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
contract=json.loads((source/'contract.json').read_text());contract.update(a=5,b_values=bs,metric='MC weighted fraction / Smin_corrected',Smin='a^2/8+9*b^2/13+a*sqrt(B)+b/2*sqrt(b^2+4*a*sqrt(B)+4*B)',normalization='MC weighted fraction / corrected Smin (Punzi Eq.8)',source_directory=str(source),source_sha256=inputs,new_fits=0,new_event_scan=False,default_changed=False);dump(out/'contract.json',contract)
books=[]
for normalized in [False,True]:
 for coord,label,axis in [('xeff_actual_percent','Weighted MC efficiency [%]','xeff'),('score_threshold','Score threshold (Prediction > threshold)','score')]:
  fig,ax=plt.subplots(figsize=(7,5))
  for b in bs:
   group=sorted([r for r in results if r['b']==b],key=lambda r:r[coord]);best=summary[str(b)];scale=best['punzi'] if normalized else 1
   curve=ax.errorbar([r[coord] for r in group],[r['punzi']/scale for r in group],yerr=[r['background_error']/scale for r in group],fmt='o-',ms=3,capsize=2,label=f"b={b:g}; max xeff={best['xeff']}%")
   color=curve[0].get_color();ax.plot(best[coord],best['punzi']/scale,'*',ms=12,color=color);ax.axvline(best[coord],ls=':',lw=1,color=color)
  ax.set(xlabel=label,ylabel='Punzi / own maximum' if normalized else 'Punzi = weighted MC fraction / Smin',title='PbPb23 v27 v1: Corrected Punzi, a=5, fine xeff grid');ax.legend(fontsize=8,loc='best');ax.text(.02,.02,'Background statistical errors only',transform=ax.transAxes,fontsize=8);fig.tight_layout()
  name=('normalized_' if normalized else '')+'xeff_scan_vs_'+axis+'.pdf';fig.savefig(out/name);plt.close(fig);books.append(out/name)
subprocess.run(['pdfunite',*map(str,books),str(out/'punzi_all.pdf')],check=True)
for n,h in inputs.items():assert sha(source/n)==h
shutil.copy2(Path(__file__),art/Path(__file__).name)
dump(out/'audit.json',dict(status='PASS',points=138,source_hashes_unchanged=True,corrected_minus_gaussian_constant_verified=True,analytic_error_vs_numerical_derivative=True,default_changed=False))
text='# v27：Corrected Punzi a=5，不同b的细扫描\n\n[四页合订本](punzi_all.pdf) · [原值vs xeff](xeff_scan_vs_xeff.pdf) · [原值vs score](xeff_scan_vs_score.pdf) · [归一化vs xeff](normalized_xeff_scan_vs_xeff.pdf) · [归一化vs score](normalized_xeff_scan_vs_score.pdf)\n\n保持count版输入，xeff15–60%每1%，46点；score横轴图是同一细网格重新排序，不是新增score均匀扫描。a=5，b=1.645、2、5，不设a=b。分母采用Punzi原文式8的经验修正式：a²/8+9b²/13+a√B+(b/2)√(b²+4a√B+4B)，公式见[contract](contract.json)。这是经验近似，不是精确Poisson检验。\n\n每张图均叠加三条曲线；原值与各自峰值归一化图分开提供，归一化仅用于比较形状。误差仅传播侧带计数，峰值归一化因子按固定常数，不表示最优点置信区间。无新拟合/事件扫描/默认变更。\n\n| b | 最优目标xeff | score cut | FOM | 99%峰高点集 |\n|---:|---:|---:|---:|---|\n'
for b in bs:
 r=summary[str(b)];text+=f"| {b:g} | {r['xeff']}% | {r['score_threshold']:.12g} | {r['punzi']:.9g} | {r['above_99percent_peak']} |\n"
text+='\nMC权重pThatreweight、保留pthat5；MC求和3.62<m<4.00，信号窗因子1；SB(3.75,3.84)/(3.90,3.99)，B=N_SB/3。共同fid18及严格Prediction>cut保持不变。修正式比完整高斯式的分母多a²/8+5b²/26，逐点已核验；即使b=5也不再与简式仅差常数比例。\n\n[数值](scan.csv) · [完整峰高集合](summary.json) · [核验](audit.json) · [Punzi原文式8](https://arxiv.org/pdf/physics/0308063)。\n'
(out/'README.md').write_text(text);print(json.dumps(summary,indent=2))
