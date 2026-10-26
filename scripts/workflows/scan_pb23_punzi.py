#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""Approved v27/v1 Punzi scans; reuse frozen xeff fits, fit 19 score cuts."""
from pathlib import Path
import argparse,json,hashlib,math,shutil,csv,re
import numpy as np
import uproot,numexpr,ROOT as R
from scipy.special import ndtr
from fit_input_contract import load_training,selection_expression
from scan_10v1_merged_punzi import SidebandFit,fingerprint,dump
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
repo=Path(__file__).resolve().parents[2];tag='X_pb23_v27_fid18_9v9_rw0_xgb_v1'
legacy=repo/'fitER/results/ml_fits'/tag/'fom_run2'
out=args.output.resolve();assert not out.exists()
R.gROOT.SetBatch(True);R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
assert R.gROOT.GetVersion()=='6.32.02'
contract=json.loads((legacy/'contract.json').read_text());oldrows=json.loads((legacy/'scan.json').read_text())
assert len(oldrows)==46 and [r['xeff'] for r in oldrows]==list(range(15,61))
assert json.loads((legacy/'independent_audit.json').read_text())['status']=='PASS'
training,oldpoints=load_training(repo,tag)
assert training['selection']==contract['selection'] and training['weight']=='pThatreweight'
protected={str(p):fingerprint(p) for p in legacy.rglob('*') if p.is_file()}
arrays={};inputs={}
# Verify raw inputs against frozen contract AND every cached column before reuse.
for key,spec in contract['inputs'].items():
 path=Path(spec['path']);fp=fingerprint(path);assert fp['sha256']==spec['sha256'],key;inputs[key]=fp
 cache=legacy/'artifacts'/f'{key}_cache.npz';cached=dict(np.load(cache))
 tree_name=training['sources'][key]['tree'] if key in training['sources'] else 'ntmix_X3872'
 branches=sorted(set(cached)|set(re.findall(r'\bB\w+\b',training['selection'])))
 with uproot.open(path) as f:raw=f[tree_name].arrays(branches,library='np')
 assert numexpr.evaluate(selection_expression(training['selection']),local_dict=raw).all()
 for k,v in cached.items():assert np.array_equal(v,raw[k]) and np.isfinite(v).all(),(key,k)
 assert np.all((cached['Prediction']>=0)&(cached['Prediction']<=1))
 if key!='data':assert np.all(cached['pThatreweight']>0)
 arrays[key]=cached
mc=arrays['mc'];data=arrays['data'];reference=arrays['reference']
D=contract['mass_range'];W=contract['window'];SB=contract['sideband_fit_range'];mean=contract['reference_mean']
assert D==[3.62,4.] and W==[3.84,3.90] and SB==[3.75,4.]
md=(mc['Bmass']>D[0])&(mc['Bmass']<D[1]);total=float(mc['pThatreweight'][md].sum())
rs=reference['Prediction'];rw=reference['pThatreweight'];order=np.argsort(rs,kind='stable');cum=np.cumsum(rw[order])/rw.sum()
for r in oldrows:
 c=float(rs[order[min(int(np.searchsorted(cum,1-r['target_efficiency'],side='left')),len(rs)-1)]])
 assert c==r['score_threshold'];assert abs(float(rw[rs>c].sum()/rw.sum())-r['achieved_efficiency'])<1e-12
 assert abs(float(mc['pThatreweight'][md&(mc['Prediction']>c)].sum())-r['mc_sumw_D'])<1e-9
out.mkdir(parents=True);art=out/'artifacts';art.mkdir()
for name in ['scan_pb23_punzi.py','plot_pb23_punzi.py','fit_input_contract.py','scan_10v1_merged_punzi.py']:
 shutil.copy2(repo/'scripts/workflows'/name,art/name)
for name in ['contract.json','scan.json','independent_audit.json']:
 shutil.copy2(legacy/name,art/('legacy_'+name))
for key,a in arrays.items():np.savez(art/(key+'_cache.npz'),**a)
dump(out/'contract.json',dict(tag=tag,metric='weighted MC fraction in D times window probability / (a/2 + sqrt(B))',a_primary=3,a_comparison=5,reference_mean=mean,reference_yield_used=False,reference_mean_origin=contract['reference_result'],mass_range=D,window=W,sideband_fit_range=SB,selection=training['selection'],weight='pThatreweight',comparison='Prediction > threshold',mc_sumw_D_no_score_cut=total,xeff_grid=list(range(15,61)),score_grid=[i/100 for i in range(5,96,5)],xeff_coordinate='achieved weighted REFERENCE MC efficiency without mass cut',xeff_fit_reuse=str(legacy),inputs=inputs,source_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in art.glob('*.py')},root_version=R.gROOT.GetVersion(),error_definition='background statistics only; no MC uncertainty or model systematic',bootstrap=0,new_data_peak_fits=0))
sbmask=(data['Bmass']>SB[0])&(data['Bmass']<SB[1])&((data['Bmass']<W[0])|(data['Bmass']>W[1]));sbfit=SidebandFit(SB,W)
root=R.TFile.Open(str(art/'score_mc_templates.root'),'RECREATE')
# MC fitting and window integration below are copied without model changes from scan_pb23_run2_fom.py.
def qual(r):return dict(status=int(r.status()),covQual=int(r.covQual()),edm=float(r.edm()),nll=float(r.minNll()))

def good(q):return q['status']==0 and q['covQual']==3 and math.isfinite(q['edm']) and q['edm']<1e-3

def template(e,mask):
 mass=R.RooRealVar('Bmass','Bmass',*D);mass.setRange('mc',*W);weight=R.RooRealVar('pThatreweight','pThatreweight',-1e6,1e6)
 ds=R.RooDataSet.from_numpy({'Bmass':mc['Bmass'][mask],'pThatreweight':mc['pThatreweight'][mask]},[mass,weight],name='mc_data',weight_name='pThatreweight')
 mu=R.RooRealVar('mean','',3.87169,3.86169,3.88169);s1=R.RooRealVar('sigma1','',.01,.001,.1);s2=R.RooRealVar('sigma2','',.005,.001,.1);frac=R.RooRealVar('fraction','',.5,.01,1.)
 g1=R.RooGaussian('g1','',mass,mu,s1);g2=R.RooGaussian('g2','',mass,mu,s2);pdf=R.RooAddPdf('pdf','',R.RooArgList(g1,g2),R.RooArgList(frac))
 opts=[R.RooFit.Save(),R.RooFit.Range('mc'),R.RooFit.Strategy(2),R.RooFit.Hesse(True),R.RooFit.SumW2Error(False),R.RooFit.PrintLevel(-1),R.RooFit.Warnings(False),R.RooFit.Verbose(False)]
 fit=pdf.fitTo(ds,*opts);attempts=[qual(fit)];root.cd();fit.Write(f'fit_xeff{e}_first')
 if not good(attempts[-1]) and math.isfinite(fit.minNll()):
  pdf.getParameters(ds).assignValueOnly(fit.floatParsFinal());fit=pdf.fitTo(ds,*opts);attempts.append(qual(fit))
 root.cd();fit.Write(f'fit_xeff{e}');ws=R.RooWorkspace(f'ws_xeff{e}');getattr(ws,'import')(pdf);ws.Write()
 nll=pdf.createNLL(ds,R.RooFit.Range('mc'));assert abs(nll.getVal()-fit.minNll())<1e-5
 result=dict(mean=mu.getVal(),sigma1=s1.getVal(),sigma2=s2.getVal(),fraction=frac.getVal(),quality_pass=good(attempts[-1]),attempts=attempts,parameter_boundary=any(min(abs(p.getVal()-p.getMin()),abs(p.getVal()-p.getMax()))<1e-4*(p.getMax()-p.getMin()) for p in [mu,s1,s2,frac]))
 return result

def acceptance(shape,limits):
 def integral(s):return (ndtr((limits[1]-mean)/s)-ndtr((limits[0]-mean)/s))/(ndtr((D[1]-mean)/s)-ndtr((D[0]-mean)/s))
 f=shape['fraction'];return float(f*integral(shape['sigma1'])+(1-f)*integral(shape['sigma2']))


def decorate(r):
 c=r['score_threshold'];sw=float(mc['pThatreweight'][md&(mc['Prediction']>c)].sum())
 assert abs(sw-r['mc_sumw_D'])<1e-9
 A=acceptance(r['mc_shape'],W);assert abs(A-r['window_acceptance'])<1e-12
 assert 0<=A<=1
 r['xeff_actual_percent']=100*float(rw[rs>c].sum()/rw.sum())
 r['mc_fraction_D']=sw/total;r['signal_window_fraction']=sw/total*A
 B=r['background'];err=r['background_error'];assert B>0
 for a in [3,5]:
  den=a/2+math.sqrt(B);r[f'punzi_a{a}']=r['signal_window_fraction']/den
  r[f'punzi_a{a}_background_error']=r['signal_window_fraction']*err/(2*math.sqrt(B)*den**2)
 return r
xrows=[]
for r in oldrows:
 q={k:r[k] for k in ['xeff','target_efficiency','score_threshold','mc_sumw_D','mc_shape','window_acceptance','background','background_error','sideband','status']}
 q['origin']='reused verified fom_run2';xrows.append(decorate(q))
dump(out/'xeff_scan.json',xrows)
srows=[]
for i in range(5,96,5):
 c=i/100;mask=md&(mc['Prediction']>c);shape=template(f'score{i:02}',mask)
 b=sbfit.fit(data['Bmass'][sbmask&(data['Prediction']>c)])
 r=dict(score_threshold=c,mc_sumw_D=float(mc['pThatreweight'][mask].sum()),mc_shape=shape,window_acceptance=acceptance(shape,W),background=b['background'],background_error=b['background_error'],sideband=b,status='PASS' if shape['quality_pass'] and b['status']=='PASS' else 'FAIL',origin='new same-model fit')
 srows.append(decorate(r));dump(out/'score_scan.json',srows)
 print('SCORE',c,'XEFF',r['xeff_actual_percent'],'PUNZI3',r['punzi_a3'],'STATUS',r['status'],flush=True)
root.Close()
summary={}
for grid,rows in [('xeff',xrows),('score',srows)]:
 valid=[r for r in rows if r['status']=='PASS'];assert valid
 summary[grid]={'points':len(rows),'valid':len(valid)}
 for a in [3,5]:
  key=f'punzi_a{a}';best=max(valid,key=lambda r:r[key])
  summary[grid][key]={'score_threshold':best['score_threshold'],'xeff_actual_percent':best['xeff_actual_percent'],'target_xeff_percent':best.get('xeff'),'value':best[key],'at_grid_boundary':best is rows[0] or best is rows[-1],'above_95percent_peak':[{'score':r['score_threshold'],'xeff_actual_percent':r['xeff_actual_percent']} for r in valid if r[key]>=.95*best[key]],'above_99percent_peak':[{'score':r['score_threshold'],'xeff_actual_percent':r['xeff_actual_percent']} for r in valid if r[key]>=.99*best[key]]}
 fields=['score_threshold','xeff_actual_percent','mc_sumw_D','mc_fraction_D','window_acceptance','signal_window_fraction','background','background_error','punzi_a3','punzi_a3_background_error','punzi_a5','punzi_a5_background_error','status']
 with (out/f'{grid}_scan.csv').open('w') as f:
  writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows({k:r[k] for k in fields} for r in rows)
dump(out/'summary.json',summary)
from plot_pb23_punzi import plot
plot(out)
for path,fp in protected.items():assert fingerprint(path)==fp,path
for key,fp in inputs.items():assert fingerprint(fp['path'])==fp,key
dump(out/'audit.json',dict(status='PASS' if all(r['status']=='PASS' for r in xrows+srows) else 'WARN_FIT_QUALITY',old_results_unchanged=True,raw_input_hash_verified=True,raw_cache_all_columns_exact=True,all_46_thresholds_exact=True,mc_sumw_closure=True,points=65,reused_fits=46,new_fits=19,protected_hashes=protected))
print('DONE',json.dumps(summary),flush=True)
