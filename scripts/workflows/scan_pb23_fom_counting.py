#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""Authorized count FOM study using frozen MC sums and reference DATA yields."""
from pathlib import Path
import argparse,json,hashlib,math,csv,shutil,subprocess,re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,o):p.write_text(json.dumps(o,indent=2,allow_nan=False)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--tag',required=True);args=ap.parse_args()
 repo=Path(__file__).resolve().parents[2];base=repo/'fitER/results/ml_fits'/args.tag;src=base/'punzi_counting';out=base/'fom_counting';assert not out.exists()
 contract=json.loads((src/'contract.json').read_text());rows=json.loads((src/'xeff_scan.json').read_text());assert [r['xeff'] for r in rows]==list(range(15,61));assert contract['mass_range']==[3.62,4.] and contract['window']==[3.84,3.90]
 assert all(r['status']=='PASS' and r['window_acceptance']==1 and abs(r['background']-(r['left_entries']+r['right_entries'])/3)<1e-10 for r in rows)
 inputs={str(p):sha(p) for p in [src/'contract.json',src/'xeff_scan.json',src/'audit.json']};refs={};lookup={r['xeff']:r for r in rows}
 for ref in [40,50]:
  p=base/f'two_peak_mc_width/xeff{ref}/fit_result.json';m=base/f'two_peak_mc_width/xeff{ref}/manifest.json';fit=json.loads(p.read_text());manifest=json.loads(m.read_text());assert fit['audit']['quality_pass'];assert manifest['threshold']==lookup[ref]['score_threshold'];assert manifest['selection']==contract['selection'];inputs[str(p)]=sha(p);inputs[str(m)]=sha(m)
  y=fit['signals']['x']['yield'];assert y['value']>0 and abs(fit['mc']['x']['sumw']/lookup[ref]['mc_sumw_D']-1)<1e-10
  refs[str(ref)]={'yield':y['value'],'yield_error':y['error'],'mc_sumw_D':lookup[ref]['mc_sumw_D'],'source':str(p),'mean_at_boundary':fit['signals']['x']['mean']['at_boundary']}
 out.mkdir();art=out/'artifacts';art.mkdir();shutil.copy2(__file__,art/Path(__file__).name)
 results=[];summaries={}
 for ref in [40,50]:
  y=refs[str(ref)]['yield'];swref=refs[str(ref)]['mc_sumw_D'];points=[]
  for r in rows:
   S=y*r['mc_sumw_D']/swref;B=r['background'];F=S/math.sqrt(S+B);err=S*r['background_error']/(2*(S+B)**1.5)
   assert abs(F*F-S*S/(S+B))<1e-12
   step=B*1e-5;g=lambda x:S/math.sqrt(S+x);assert abs(abs((g(B+step)-g(B-step))/(2*step))*r['background_error']-err)<1e-9
   q=dict(reference_xeff=ref,xeff=r['xeff'],xeff_actual_percent=r['xeff_actual_percent'],score_threshold=r['score_threshold'],mc_sumw_D=r['mc_sumw_D'],S=S,B=B,FOM=F,background_error_only=err,punzi_a5=r['punzi_a5']);points.append(q)
  assert abs(next(q for q in points if q['xeff']==ref)['S']-y)<1e-10
  best=max(points,key=lambda r:r['FOM']);summaries[str(ref)]={**best,'above_99percent_peak':[q['xeff'] for q in points if q['FOM']>=.99*best['FOM']],'above_95percent_peak':[q['xeff'] for q in points if q['FOM']>=.95*best['FOM']]};results+=points
 dump(out/'scan.json',results);dump(out/'summary.json',summaries)
 with (out/'scan.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
 label=args.tag.replace('X_pb23_','PbPb23 ').replace('_rw0_xgb_v1','')
 for coord,xlabel in [('xeff_actual_percent','Weighted MC efficiency [%]'),('score_threshold','Score threshold (Prediction > threshold)')]:
  suffix='xeff' if coord=='xeff_actual_percent' else 'score'
  for normalized in [False,True]:
   fig,ax=plt.subplots(figsize=(7,5))
   for ref in [40,50]:
    pts=sorted([r for r in results if r['reference_xeff']==ref],key=lambda r:r[coord]);best=summaries[str(ref)];norm=best['FOM'] if normalized else 1
    line=ax.errorbar([r[coord] for r in pts],[r['FOM']/norm for r in pts],yerr=None if normalized else [r['background_error_only'] for r in pts],fmt='o-',ms=3,capsize=2,label=f"FOM ref {ref}%: max {best['xeff']}%")
    ax.plot(best[coord],best['FOM']/norm,'*',ms=12,color=line[0].get_color())
   if normalized:
    pts=sorted(rows,key=lambda r:r[coord]);peak=max(r['punzi_a5'] for r in pts);ax.plot([r[coord] for r in pts],[r['punzi_a5']/peak for r in pts],'--',label='Punzi a=5 / its maximum')
   ax.set(xlabel=xlabel,ylabel='Metric / own maximum' if normalized else 'S / sqrt(S+B)',title=label);ax.legend(fontsize=8,loc='best');fig.tight_layout();fig.savefig(out/f'{"comparison" if normalized else "fom"}_vs_{suffix}.pdf');plt.close(fig)
 pdfs=[out/f'{kind}_vs_{coord}.pdf' for kind in ['fom','comparison'] for coord in ['xeff','score']];subprocess.run(['pdfunite',*map(str,pdfs),str(out/'fom_all.pdf')],check=True)
 assert int(re.search(r'^Pages:\s+(\d+)',subprocess.check_output(['pdfinfo',str(out/'fom_all.pdf')],text=True),re.M)[1])==4
 for p,h in inputs.items():assert sha(Path(p))==h
 dump(out/'contract.json',dict(tag=args.tag,metric='S/sqrt(S+B)',S='reference first-stage DATA X yield * MC_sumw_D(c)/MC_sumw_D(reference)',B='(N(3.75<m<3.84)+N(3.90<m<3.99))/3',signal_window_factor=1,mass_range=[3.62,4.],selection=contract['selection'],mc_apply_selection=contract['mc_apply_selection'],references=refs,inputs=inputs,grid=list(range(15,61)),errors='plotted errors: background Poisson only; reference yield uncertainty common/correlated, recorded but not propagated; no MC or model uncertainties',new_fits=0,bootstrap=0))
 dump(out/'audit.json',dict(status='PASS',points=92,reference_yield_closure=True,reference_threshold_selection_mc_sumw_match=True,sideband_width_ratio=True,independent_squared_formula=True,error_derivative=True,inputs_unchanged=True,pdf_pages=4))
 (out/'README.md').write_text('# '+args.tag+': counting FOM\n\n[四页图册](fom_all.pdf) · [数值CSV](scan.csv) · [契约](contract.json) · [验收](audit.json)\n\nS = 参考点第一步DATA X yield × MC加权通过量/参考点MC加权通过量；MC窗口因子1。B = 对称侧带(3.75,3.84)/(3.90,3.99)事例数/3。FOM=S/sqrt(S+B)。分别以40%、50%为参考，xeff15–60%每1%；两种横轴是同一46点。\n\n'+''.join(f"- 参考{ref}%：最大xeff={summaries[str(ref)]['xeff']}%，FOM={summaries[str(ref)]['FOM']:.6f}；99%峰高点{summaries[str(ref)]['above_99percent_peak']}。\n" for ref in [40,50])+'\n误差条只含背景计数统计；参考yield误差已记录，其相关传播未计入。归一化比较图每条曲线除以自身最大值，不比较不同FOM的绝对数值。峰高点集不是置信区间，FOM不是拟合显著性。无新拟合、bootstrap；全部复用冻结计数，输入hash不变。\n')
 dump(out/'delivery_hashes.json',{str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()});print(args.tag,json.dumps(summaries),flush=True)
if __name__=='__main__':main()
