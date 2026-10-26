#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""User-approved symmetric-sideband counting Punzi, with signal window factor 1."""
from pathlib import Path
import argparse,json,hashlib,math,csv,shutil
import numpy as np
from plot_pb23_punzi_a5 import plot
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
repo=Path(__file__).resolve().parents[2];tag='X_pb23_v27_fid18_9v9_rw0_xgb_v1';source=repo/'fitER/results/ml_fits'/tag/'punzi_fom';out=args.output.resolve();assert not out.exists()
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
protected={str(p):sha(p) for p in source.rglob('*') if p.is_file()}
for name,h in json.loads((source/'delivery_hashes.json').read_text()).items():assert sha(source/name)==h
previous=json.loads((source/'contract.json').read_text())
for spec in previous['inputs'].values():assert sha(spec['path'])==spec['sha256']
a={k:dict(np.load(source/'artifacts'/f'{k}_cache.npz')) for k in ['data','mc','reference']}
mc=a['mc'];data=a['data'];reference=a['reference'];md=(mc['Bmass']>3.62)&(mc['Bmass']<4.0);total=float(mc['pThatreweight'][md].sum());assert abs(total-previous['mc_sumw_D_no_score_cut'])<1e-10
# Old fid18 sample retains all its delivered pThat weight categories.
# Reuse the verified preselected fid18 cache. Only Bmass windows and score masks change here.
left=(data['Bmass']>3.75)&(data['Bmass']<3.84);right=(data['Bmass']>3.90)&(data['Bmass']<3.99)
ratio=1/3;out.mkdir(parents=True);art=out/'artifacts';art.mkdir()
for p in [Path(__file__),repo/'scripts/workflows/plot_pb23_punzi_a5.py']:shutil.copy2(p,art/p.name)
contract={**previous,'plot_label':'PbPb23 v27 9v9 v1: counting','metric':'MC weighted fraction in D / (a/2 + sqrt((N_left+N_right)/3))','a_primary':5,'signal_window_factor':1.,'reference_mean_used':False,'reference_yield_used':False,'sideband_windows':[[3.75,3.84],[3.90,3.99]],'sideband_comparison':'strict mass boundaries; Prediction > threshold','background_width_ratio':ratio,'background_method':'symmetric sideband counts, no fit','cache_source':str(source/'artifacts'),'source_input_contract_sha256':sha(source/'contract.json'),'cache_sha256':{k:sha(source/'artifacts'/f'{k}_cache.npz') for k in a},'error_definition':'sideband Poisson count only; no MC statistics or model systematic','source_hashes':{p.name:sha(p) for p in art.glob('*.py')}}
for key in ['reference_mean','reference_mean_origin','reference_xeff','sideband_fit_range','a_comparison']:contract.pop(key,None)
dump(out/'contract.json',contract);summary={};allrows=[]
for grid in ['xeff','score']:
 old=json.loads((source/f'{grid}_scan.json').read_text());rows=[]
 for r in old:
  c=r['score_threshold'];score=data['Prediction']>c;nl=int(np.count_nonzero(left&score));nr=int(np.count_nonzero(right&score));n=nl+nr;assert n>0
  B=n/3;err=math.sqrt(n)/3;sw=float(mc['pThatreweight'][md&(mc['Prediction']>c)].sum());fraction=sw/total
  assert abs(fraction-r['mc_fraction_D'])<1e-12
  eff=100*float(reference['pThatreweight'][reference['Prediction']>c].sum()/reference['pThatreweight'].sum());assert abs(eff-r['xeff_actual_percent'])<1e-10
  row=dict(score_threshold=c,xeff_actual_percent=eff,xeff=r.get('xeff'),mc_sumw_D=sw,mc_fraction_D=fraction,window_acceptance=1.,signal_window_fraction=fraction,left_entries=nl,right_entries=nr,sideband_entries=n,background=B,background_error=err,previous_background=r['background'],previous_window_acceptance=r['window_acceptance'],status='PASS')
  for av in [5]:
   den=av/2+math.sqrt(B);row[f'punzi_a{av}']=fraction/den;row[f'punzi_a{av}_background_error']=fraction*err/(2*math.sqrt(B)*den**2);row[f'previous_punzi_a{av}']=r[f'punzi_a{av}']
  rows.append(row)
 dump(out/f'{grid}_scan.json',rows)
 with (out/f'{grid}_scan.csv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 summary[grid]={'points':len(rows)}
 for av in [5]:
  key=f'punzi_a{av}';best=max(rows,key=lambda r:r[key]);summary[grid][key]={**{k:best[k] for k in ['score_threshold','xeff_actual_percent','xeff']},'value':best[key],'at_grid_boundary':best is rows[0] or best is rows[-1],'above_95percent_peak':[r['xeff'] if grid=='xeff' else r['score_threshold'] for r in rows if r[key]>=.95*best[key]],'above_99percent_peak':[r['xeff'] if grid=='xeff' else r['score_threshold'] for r in rows if r[key]>=.99*best[key]]}
 allrows+=rows
dump(out/'summary.json',summary);plot(out)
# Independent checks use sorted interval counts (rather than mass/score boolean intersections),
# fsum instead of numpy sums, and a numerical derivative for uncertainty propagation.
checks=[]
for r in allrows:
 c=r['score_threshold'];m=np.sort(data['Bmass'][data['Prediction']>c])
 counts=[int(np.searchsorted(m,h,side='left')-np.searchsorted(m,l,side='right')) for l,h in [(3.75,3.84),(3.90,3.99)]]
 assert counts==[r['left_entries'],r['right_entries']]
 fraction=math.fsum(mc['pThatreweight'][md&(mc['Prediction']>c)])/math.fsum(mc['pThatreweight'][md]);b=sum(counts)*(3.90-3.84)/((3.84-3.75)+(3.99-3.90));assert abs(b-r['background'])<1e-8
 for av in [5]:
  fun=lambda bb:fraction/(av/2+math.sqrt(bb));assert abs(fun(b)-r[f'punzi_a{av}'])<1e-12
  step=b*1e-5;error=abs((fun(b+step)-fun(b-step))/(2*step))*math.sqrt(sum(counts))/3;assert abs(error-r[f'punzi_a{av}_background_error'])<1e-10
 checks.append(c)
assert len(checks)==65
for name,h in protected.items():assert sha(name)==h
for spec in previous['inputs'].values():assert sha(spec['path'])==spec['sha256']
dump(out/'audit.json',dict(status='PASS',points=65,independent_sorted_mass_counts=True,independent_fsum_and_width_ratio=True,numerical_derivative_error_verified=True,original_results_and_inputs_unchanged=True,signal_window_factor_exactly_one=True,new_fits=0,protected_hashes=protected))
print(json.dumps(summary,indent=2))
