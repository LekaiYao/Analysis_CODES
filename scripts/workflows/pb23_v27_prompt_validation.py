#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""v27/v1 xeff50 signed BLxy < 0.01 cm check: refit MC/DATA and sPlot."""
from pathlib import Path
import argparse, json, hashlib, shutil, math
import numpy as np
import uproot, ROOT as R
from ppref_psi2s_pthat_validation import fingerprint,dump,stats,record
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',required=True,type=Path);args=ap.parse_args()
repo=Path(__file__).resolve().parents[2];source=repo/'plotER/validation_archive/results/v27_xeff50_splot';out=args.output.resolve();assert not out.exists()
R.gROOT.SetBatch(True);assert R.gROOT.GetVersion()=='6.32.02';R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
parent=json.loads((source/'summary.json').read_text());variables=parent['variables'];cut=parent['score_threshold'];fixed_width=True
out.mkdir(parents=True);art=out/'artifacts';art.mkdir();protected={'parent_summary':source/'summary.json','script':Path(__file__)}
arrays={};records={};counts={};handles=[];keep=[]
for key,tree,name in [('data','ntmix_sweights','data_sweights.root'),('x','ntmix_X3872','x_mc.root'),('psi','ntmix_PSI2S','psi_mc.root')]:
 path=source/'artifacts'/name;protected[key]=path
 with uproot.open(path) as f:a=f[tree].arrays(library='np')
 mask=np.isfinite(a['BLxy'])&(a['BLxy']<.01)
 arrays[key]={n:v[mask] for n,v in a.items() if not n.endswith('_sWeight')}
 counts[key]={'inclusive':len(mask),'retained':int(mask.sum()),'rejected_nonfinite':int((~np.isfinite(a['BLxy'])).sum())}
 if key!='data':counts[key]['pthat_weighted_efficiency']=float(a['pThatreweight'][mask].astype(float).sum()/a['pThatreweight'].astype(float).sum())
 b=arrays[key]; assert np.all(b['Prediction']>cut)
 fitmask=(b['Bmass']>3.62)&(b['Bmass']<4.0);cache=art/(key+'_fit_cache.root');columns=['Bmass']+(['pThatreweight'] if key!='data' else [])
 with uproot.recreate(cache) as f:f.mktree(tree,{n:b[n].dtype for n in columns});f[tree].extend({n:b[n][fitmask] for n in columns})
 records[key]={'cache':str(cache),'tree':tree,'entries':int(fitmask.sum())}
provenance={k:fingerprint(v) for k,v in protected.items()};dump(out/'provenance.json',provenance)
shutil.copy2(__file__,art/Path(__file__).name)
mass=R.RooRealVar('Bmass','Bmass',3.62,4.0);mass.setRange('all',3.62,4.0);weight=R.RooRealVar('pThatreweight','pThatreweight',-1e6,1e6)
def ds(key):
 info=records[key];f=R.TFile.Open(info['cache']);handles.append(f);args=R.RooArgSet(mass)
 if key!='data':args.add(weight)
 opts=[R.RooFit.Import(f.Get(info['tree']))]+([R.RooFit.WeightVar('pThatreweight')] if key!='data' else []);d=R.RooDataSet(key+'_data',key+'_data',args,*opts);keep.append(d);assert d.numEntries()==info['entries'];return d
data=ds('data');mc={k:ds(k) for k in ['x','psi']};params={};pdfs={};mcresults={};fitresults={}
def var(n,v,lo,hi):
 z=R.RooRealVar(n,n,v,lo,hi);keep.append(z);return z
retry_history=[]
def fit_with_one_restart(pdf,d,opts):
 def snapshot(q):
  return {'status':q.status(),'covQual':q.covQual(),'edm':q.edm(),'nll':q.minNll(),'history':[(q.statusLabelHistory(i),q.statusCodeHistory(i)) for i in range(q.numStatusHistory())],'parameters':{v.GetName():{'value':v.getVal(),'error':v.getError()} for v in q.floatParsFinal()}}
 q=pdf.fitTo(d,*opts);attempts=[snapshot(q)]
 if math.isfinite(q.minNll()) and (q.status()!=0 or q.covQual()!=3 or not math.isfinite(q.edm()) or q.edm()>=1e-3):
  pars=pdf.getParameters(d);pars.assignValueOnly(q.floatParsFinal())
  q=pdf.fitTo(d,*opts);attempts.append(snapshot(q))
 retry_history.append({'pdf':pdf.GetName(),'attempts':attempts})
 (out/'retry_history.json').write_text(json.dumps(retry_history,indent=2)+'\n')
 return q
def fit(pdf,d,window,weighted=False):
 opts=[R.RooFit.Save(),R.RooFit.Range(window),R.RooFit.Strategy(2),R.RooFit.Hesse(True),R.RooFit.PrintLevel(-1),R.RooFit.Warnings(False),R.RooFit.Verbose(False),R.RooFit.SumW2Error(False) if weighted else R.RooFit.Extended(True)];r=fit_with_one_restart(pdf,d,opts);keep.append(r);return r
for key,nom,lo,hi in [('x',3.87169,3.84,3.9),('psi',3.68610,3.65610,3.71610)]:
 mu=var(key+'_mc_mean',nom,nom-.01,nom+.01);s1=var(key+'_sigma1',.01,.001,.1);s2=var(key+'_sigma2',.005,.001,.1);frac=var(key+'_fraction',.5,.01,1.)
 g1=R.RooGaussian(key+'_mc_g1','',mass,mu,s1);g2=R.RooGaussian(key+'_mc_g2','',mass,mu,s2);pdf=R.RooAddPdf(key+'_mc_pdf','',R.RooArgList(g1,g2),R.RooArgList(frac));keep.extend([g1,g2,pdf]);mass.setRange(key+'_mc_range',lo,hi);r=fit(pdf,mc[key],key+'_mc_range',True);fitresults[key+'_mc']=r
 mcresults[key]={'mean':mu.getVal(),'sigma1':s1.getVal(),'sigma2':s2.getVal(),'fraction':frac.getVal(),'range':[lo,hi],'sumw':mc[key].sumEntries(),'sumw2':float((arrays[key]['pThatreweight'][(arrays[key]['Bmass']>3.62)&(arrays[key]['Bmass']<4.0)]**2).sum())}
 for p in [mu,s1,s2,frac]:p.setConstant(True)
 mean=var(key+'_mean',mu.getVal(),nom-.01,nom+.01);scale=var(key+'_scale',1.,.5,1.5);ss1=R.RooProduct(key+'_scaled_sigma1','',R.RooArgList(scale,s1));ss2=R.RooProduct(key+'_scaled_sigma2','',R.RooArgList(scale,s2));dg1=R.RooGaussian(key+'_g1','',mass,mean,ss1);dg2=R.RooGaussian(key+'_g2','',mass,mean,ss2);sig=R.RooAddPdf(key+'_signal','',R.RooArgList(dg1,dg2),R.RooArgList(frac));near=int(np.sum(abs(arrays['data']['Bmass']-nom)<.005));yieldv=var(key+'_yield',max(0.,.4*near),0.,max(10.,2.*near));params[key]={'mean':mean,'scale':scale,'yield':yieldv};pdfs[key]=sig;keep.extend([ss1,ss2,dg1,dg2,sig]);pdfs[key+'_mc']=pdf
if fixed_width:
 for pp in params.values():pp['scale'].setVal(1.);pp['scale'].setConstant(True)
n=data.numEntries();a0=var('a0',-.35,-2,2);a1=var('a1',-.05,-2,2);nb=var('nbkg',.7*n,.1*n,n);bg=R.RooChebychev('backgroundPdf','',mass,R.RooArgList(a0,a1));model=R.RooAddPdf('model','',R.RooArgList(pdfs['x'],pdfs['psi'],bg),R.RooArgList(params['x']['yield'],params['psi']['yield'],nb));keep.extend([bg,model]);alt=fit(model,data,'all');fitresults['alt']=alt;allpars=model.getParameters(data);snapshot=allpars.snapshot()
def quality(r):return {'status':r.status(),'covQual':r.covQual(),'edm':r.edm(),'nll':r.minNll()}
def boundary(p):return min(abs(p.getVal()-p.getMin()),abs(p.getVal()-p.getMax()))<=1e-4*(p.getMax()-p.getMin())

# Numerical gate for MC templates and full DATA fit; at most the existing one restart.
for label,r in fitresults.items():
 assert r.status()==0 and r.covQual()==3 and math.isfinite(r.edm()) and r.edm()<1e-3,(label,record(r))
 assert all(not q['near_boundary'] for q in record(r)['parameters'].values()),(label,record(r))
pts=[-1.,1.];aa=a0.getVal();bb=a1.getVal()
if bb and -1 < -aa/(4*bb) < 1:pts.append(-aa/(4*bb))
assert min(1+aa*z+bb*(2*z*z-1) for z in pts)>0
full_fit=record(alt);dump(out/'fit_result.json',{'data':full_fit,'mc':{k:record(fitresults[k+'_mc']) for k in ['x','psi']},'mc_shapes':mcresults,'shape_strategy':'mc_width_mean_float','extra_selection':'BLxy < 0.01 cm (signed)'})
w=R.RooWorkspace('ws_nominal');imp=getattr(w,'import');imp(data);imp(model)
fout=R.TFile(str(out/'fit_workspace.root'),'RECREATE');w.Write();alt.Write('fit_result_alt');fout.Close()
# Accepted mass-fit template, with selection annotation replacing unrequested significance.
plot=(repo/'fitER/models/TwoPeakPlots.C').read_text();needle='stats.AddText(Form("Z_{PL,X}=%.3f", z));';assert plot.count(needle)==1
plot=plot.replace(needle,'stats.AddText("L_{xy} < 0.01 cm");')
(art/'TwoPeakPlots.C').write_text(plot);R.gSystem.SetBuildDir(str(art),True);assert R.gSystem.CompileMacro(str(art/'TwoPeakPlots.C'),'k')>0
R.PlotTwoData(str(out/'mass_fit.pdf'),data,model,pdfs['x'],pdfs['psi'],bg,mass,alt,0.)
for key,label in [('x','X'),('psi','#psi(2S)')]:
 info=mcresults[key];R.PlotPeakMC(str(out/(key+'_mc_template_fit.pdf')),mc[key],pdfs[key+'_mc'],mass,fitresults[key+'_mc'],info['mean'],info['sigma1'],info['sigma2'],info['fraction'],*info['range'],label)
model=w.pdf('model');data=w.data('data_data');mass=w.var('Bmass');pars=model.getParameters(data);d=arrays['data']
ykeys=['x_yield','psi_yield','nbkg'];yields=[w.var(k) for k in ykeys]
for p in pars:p.setConstant(p.GetName() not in ykeys)
fixed={p.GetName():p.getVal() for p in pars if p.isConstant()};assert fixed['x_scale']==fixed['psi_scale']==1.
opts=[R.RooFit.Save(),R.RooFit.Extended(True),R.RooFit.Range('all'),R.RooFit.Strategy(2),R.RooFit.Hesse(True),R.RooFit.PrintLevel(-1)]
good=lambda q:q.status()==0 and q.covQual()==3 and math.isfinite(q.edm()) and q.edm()<1e-3
fit=model.fitTo(data,*opts);history=[record(fit)]
if not good(fit) and math.isfinite(fit.minNll()):pars.assignValueOnly(fit.floatParsFinal());fit=model.fitTo(data,*opts);history.append(record(fit))
dump(out/'splot_retry_history.json',history);assert good(fit),history;assert set(p.GetName() for p in fit.floatParsFinal())==set(ykeys)
assert all(not p['near_boundary'] for p in record(fit)['parameters'].values())
pre={p.GetName():p.getVal() for p in yields}
# ROOT 6.32 SPlot sets yields to zero BEFORE widening a positive lower bound.
# Prepare the [0,1] probing domain first, avoiding a clamped background at event 0.
original_min={p.GetName():p.getMin() for p in yields}
for p in yields:p.setMin(0.)
spl=R.RooStats.SPlot('splot','splot',data,model,R.RooArgList(*yields),R.RooArgSet(),True,False,'',R.RooFit.Strategy(2),R.RooFit.Hesse(True),R.RooFit.PrintLevel(-1),R.RooFit.Range('all'))
sw=np.array([[spl.GetSWeight(i,k) for k in ykeys] for i in range(data.numEntries())]);assert np.isfinite(sw).all()
for p in yields:
 assert p.getVal()>original_min[p.GetName()]
 p.setMin(original_min[p.GetName()])
for name,val in fixed.items():assert w.var(name).isConstant() and w.var(name).getVal()==val
assert all(abs(p.getVal()-pre[p.GetName()])<.01*fit.floatParsFinal().find(p.GetName()).getError() for p in yields)
components=[w.pdf(k) for k in ['x_signal','psi_signal','backgroundPdf']];obs=R.RooArgSet(mass);P=[]
for value in d['Bmass']:
 mass.setVal(float(value));P.append([p.getVal(obs) for p in components])
P=np.array(P);nv=np.array([p.getVal() for p in yields]);den=P@nv;assert np.all(den>0);A=P/den[:,None];information=A.T@A;V=np.linalg.inv(information);analytic=A@V.T
L=np.array([[data.get(i).getRealValue('L_'+k) for k in ykeys] for i in range(data.numEntries())]);assert np.max(abs(L-P))<1e-8
assert np.linalg.eigvalsh(information).min()>0;assert np.max(abs(sw-analytic))<1e-8
assert np.max(abs(sw.sum(axis=0)-nv))<1e-5*max(nv);assert np.max(abs(sw.sum(axis=1)-1))<1e-3
assert np.max(abs(A.sum(axis=0)-1))<1e-4
for i,name in enumerate(['x_sWeight','psi_sWeight','background_sWeight']):d[name]=sw[:,i]
for key,tree,name in [('data','ntmix_sweights','data_sweights.root'),('x','ntmix_X3872','x_mc.root'),('psi','ntmix_PSI2S','psi_mc.root')]:
 with uproot.recreate(art/name) as ff:ff.mktree(tree,{k:v.dtype for k,v in arrays[key].items()});ff[tree].extend(arrays[key])
np.savez(art/'sweight_algebra.npz',component_pdf=P,yields=nv,covariance=V,sweights=sw)
rf=R.TFile(str(art/'splot_workspace.root'),'RECREATE');rf.cd();w.Write();data.Write('splot_data');fit.Write('fit_yields');rf.Close()

results={k:stats(sw[:,i]) for i,k in enumerate(['x','psi','background'])}
summary={'status':'SPLOT_PASS','tag':parent['tag'],'xeff':50,'score_threshold':cut,'selection':'('+parent['selection']+') and (BLxy < 0.01)','extra_selection':'BLxy < 0.01 cm (signed, no abs)','data_mass_range':[3.62,4.0],'mc_validation_mass_cut':None,'mc_weight':'pThatreweight','fit_source':str(out),'fit':record(fit),'full_shape_fit':full_fit,'splot_yields':dict(zip(ykeys,nv.tolist())),'sweights':results,'mc':{k:stats(arrays[k]['pThatreweight'].astype(float)) for k in ['x','psi']},'cutflow':counts,'variables':variables,'fixed_parameters':fixed,'checks':{'analytic_sweight_max_difference':float(np.max(abs(sw-analytic))),'sumw_minus_yield':(sw.sum(axis=0)-nv).tolist()},'limitations':['Prompt-enriched, residual nonprompt not subtracted.','sPlot assumes mass independence after ML and Lxy selections.','Fixed shape uncertainty not propagated; no bootstrap.']}
dump(out/'summary.json',summary)
inclusive=json.loads((Path(parent['fit_source'])/'fit_result.json').read_text())
comparison={}
for key in ['x','psi']:
 inc=inclusive['signals'][key]['yield'];p=full_fit['parameters'][key+'_yield']
 comparison[key]={'inclusive_full_fit':inc,'prompt_enriched_full_fit':p,'retention_central':p['value']/inc['value'],'yield_change':p['value']-inc['value'],'inclusive_splot':parent['sweights'][key],'prompt_enriched_splot':results[key],'mc_pthat_efficiency':counts[key]['pthat_weighted_efficiency']}
dump(out/'yield_comparison.json',{'species':comparison,'ratio_uncertainty':'not assigned: overlapping fits with refitted shapes; individual fit uncertainties shown, no independent-sample ratio error','interpretation':'selection retention, not prompt purity or prompt fraction'})
for k,path in protected.items():assert fingerprint(path)==provenance[k]
print('PROMPT SPLOT PASS',json.dumps(comparison),flush=True)
