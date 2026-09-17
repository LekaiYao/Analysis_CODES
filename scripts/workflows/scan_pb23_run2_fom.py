#!/usr/bin/env python3
"""Approved v27 Run-2 counting FOM, frozen 50% DATA normalization; no DATA peak scan."""
from pathlib import Path
import json, hashlib, csv, math, re, shutil
import numpy as np
import uproot, numexpr, ROOT as R
from scipy.special import ndtr
from fit_input_contract import load_training, selection_expression
from scan_10v1_merged_punzi import SidebandFit, fingerprint, dump
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R.gROOT.SetBatch(True); R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
repo=Path(__file__).resolve().parents[2];tag='X_pb23_v27_fid18_9v9_rw0_xgb_v1'
base=repo/'fitER/results/ml_fits'/tag;out=base/'fom_run2';assert not out.exists()
assert R.gROOT.GetVersion()=='6.32.02'
training,oldpoints=load_training(repo,tag);assert training['weight']=='pThatreweight'
reference=base/'checks/two_peak_mc_width/xeff50';ref=json.loads((reference/'fit_result.json').read_text());rm=json.loads((reference/'manifest.json').read_text())
assert ref['audit']['quality_pass'] and rm['selection']==training['selection']
N=ref['signals']['x']['yield']['value'];Nerr=ref['signals']['x']['yield']['error'];mean=ref['signals']['x']['mean']['value'];c0=rm['threshold']
D=[3.62,4.0];W=[3.84,3.90];SB=[3.75,4.0];grid=list(range(15,61))
refpath=repo.parent/'XGBoost/output/selected'/tag/'REFERENCE_MC_with_score.root'
inputs={};arrays={}
for key,spec in {**training['sources'],'reference':{'path':str(refpath),'tree':'ntmix_X3872'}}.items():
 path=Path(spec['path']);inputs[key]=fingerprint(path)
 branches=sorted(set(['Bmass','Prediction']+re.findall(r'\bB\w+\b',training['selection'])+(['pThatreweight'] if key!='data' else [])))
 with uproot.open(path) as f:a=f[spec['tree']].arrays(branches,library='np')
 assert all(np.isfinite(v).all() for v in a.values())
 fid=numexpr.evaluate(selection_expression(training['selection']),local_dict=a);assert fid.all(),key
 arrays[key]={k:np.asarray(a[k],dtype=float) for k in ['Bmass','Prediction']+(['pThatreweight'] if key!='data' else [])}
 if key!='data':assert np.all(arrays[key]['pThatreweight']>0)
 assert fingerprint(path)==inputs[key]
rs=arrays['reference']['Prediction'];rw=arrays['reference']['pThatreweight'];order=np.argsort(rs,kind='stable');cum=np.cumsum(rw[order])/rw.sum()
def threshold(e):
 i=min(int(np.searchsorted(cum,1-e,side='left')),len(rs)-1);c=float(rs[order[i]])
 return dict(target_efficiency=e,score_threshold=c,achieved_efficiency=float(rw[rs>c].sum()/rw.sum()))
for point in oldpoints:
 new=threshold(point['target_efficiency']);assert new['score_threshold']==point['score_threshold'];assert abs(new['achieved_efficiency']-point['achieved_efficiency'])<1e-12
points=[threshold(e/100) for e in grid];assert points[grid.index(50)]['score_threshold']==c0
out.mkdir();art=out/'artifacts';art.mkdir();shutil.copy2(__file__,art)
for name in ['fit_input_contract.py','scan_10v1_merged_punzi.py']:shutil.copy2(repo/'scripts/workflows'/name,art)
shutil.copy2(repo.parent/'XGBoost/utils/score_thresholds.py',art/'upstream_score_thresholds.py')
for key,a in arrays.items():np.savez(art/(key+'_cache.npz'),**a)
contract=dict(tag=tag,grid=grid,selection=training['selection'],reference_yield=N,reference_error=Nerr,reference_mean=mean,reference_cut=c0,reference_result=fingerprint(reference/'fit_result.json'),reference_manifest=fingerprint(reference/'manifest.json'),inputs=inputs,training=training,mass_range=D,window=W,sideband_fit_range=SB,background='conditional unbinned C2, [-2,2], positive over full range',mc_shape='same-cut weighted unbinned double Gaussian, [3.84,3.90], Strategy2 SumW2Error(false); one permitted warm restart',window_shape_mean='frozen 50% DATA mean; MC widths/fraction, scale=1',normalization='MC efficiency in D relative to c0; window PDF integral normalized in D',interpretation='inclusive DATA reference transported with prompt MC efficiency; uncalibrated counting proxy',bootstrap=0,root_version=R.gROOT.GetVersion(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
dump(out/'contract.json',contract)
supp=dict(schema_version=1,train_tag=tag,parent_thresholds_sha256=training['thresholds_sha256'],selection=training['selection'],weight_branch=['pThatreweight'],comparison='Prediction > threshold',existing_thresholds_exactly_reproduced=len(oldpoints),reference_file=str(refpath),reference_sha256=inputs['reference']['sha256'],thresholds=points)
dump(out/'thresholds.json',supp)
mc=arrays['mc'];md=(mc['Bmass']>D[0])&(mc['Bmass']<D[1]);norm=float(mc['pThatreweight'][md&(mc['Prediction']>c0)].sum());total=float(mc['pThatreweight'][md].sum())
assert np.isclose(norm,ref['mc']['x']['sumw'],rtol=1e-12)
data=arrays['data'];sbmask=(data['Bmass']>SB[0])&(data['Bmass']<SB[1])&((data['Bmass']<W[0])|(data['Bmass']>W[1]));sbfit=SidebandFit(SB,W)
root=R.TFile.Open(str(art/'mc_templates.root'),'RECREATE');keep=[]
def qual(r):return dict(status=int(r.status()),covQual=int(r.covQual()),edm=float(r.edm()),nll=float(r.minNll()))
def good(q):return q['status']==0 and q['covQual']==3 and math.isfinite(q['edm']) and q['edm']<1e-3
def template(e,mask):
 if e==50:return dict(ref['mc']['x'],quality_pass=True,reused_reference=True,attempts=[])
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
rows=[]
for e,p in zip(grid,points):
 c=p['score_threshold'];mask=md&(mc['Prediction']>c);sw=float(mc['pThatreweight'][mask].sum());shape=template(e,mask);A=acceptance(shape,W);b=sbfit.fit(data['Bmass'][sbmask&(data['Prediction']>c)])
 signal=N*sw/norm*A;B=b['background'];fom=signal/math.sqrt(signal+B)
 r=dict(xeff=e,**p,mc_efficiency_D=sw/total,relative_mc_efficiency=sw/norm,mc_sumw_D=sw,mc_shape=shape,window_acceptance=A,signal=signal,background=B,background_error=b['background_error'],fom=fom,fom_background_error=signal/(2*(signal+B)**1.5)*b['background_error'],fom_low_reference=(signal*(N-Nerr)/N)/math.sqrt(signal*(N-Nerr)/N+B),fom_high_reference=(signal*(N+Nerr)/N)/math.sqrt(signal*(N+Nerr)/N+B),sideband=b,signal_tail_left=N*sw/norm*acceptance(shape,[SB[0],W[0]]),signal_tail_right=N*sw/norm*acceptance(shape,[W[1],SB[1]]),status='PASS' if shape['quality_pass'] and b['status']=='PASS' else 'FAIL')
 rows.append(r);dump(out/'scan.json',rows);print('POINT',e,'FOM',fom,'S',signal,'B',B,r['status'],flush=True)
root.Close()
valid=[r for r in rows if r['status']=='PASS'];assert valid
best=max(valid,key=lambda r:r['fom']);plateau=[r['xeff'] for r in valid if r['fom']>=.95*best['fom']]
selected=dict(xeff=best['xeff'],score_threshold=best['score_threshold'],fom=best['fom'],signal=best['signal'],background=best['background'],peak_at_scan_boundary=best['xeff'] in [grid[0],grid[-1]],points_above_95percent_peak=plateau,low_reference_peak=max(valid,key=lambda r:r['fom_low_reference'])['xeff'],high_reference_peak=max(valid,key=lambda r:r['fom_high_reference'])['xeff'],selection='largest valid Run-2 counting FOM; frozen before new selected-point DATA fit',scan_sha256=hashlib.sha256((out/'scan.json').read_bytes()).hexdigest())
dump(out/'selected_point.json',selected)
fields=['xeff','score_threshold','achieved_efficiency','mc_efficiency_D','relative_mc_efficiency','window_acceptance','signal','background','background_error','fom','fom_background_error','fom_low_reference','fom_high_reference','status']
with (out/'scan.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows({k:r[k] for k in fields} for r in rows)
from plot_pb23_run2_fom import plot
plot(out)
for key,spec in inputs.items():assert fingerprint(spec['path'])==spec
assert all(points[i]['score_threshold']>=points[i+1]['score_threshold'] for i in range(len(points)-1))
dump(out/'audit.json',dict(status='PASS' if len(valid)==len(rows) else 'WARN_FIT_QUALITY',points=len(rows),valid_points=len(valid),input_hashes_unchanged=True,existing_threshold_closure=len(oldpoints),reference_mc_sumw_closure=True,reference_signal_closure=abs(rows[grid.index(50)]['signal']-N*rows[grid.index(50)]['window_acceptance'])<1e-10,selected_point_frozen=True))
print('SELECTED',json.dumps(selected),flush=True)
