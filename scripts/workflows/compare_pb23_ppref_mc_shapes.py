"""PbPb23 prompt MC with the ppRef ordinary precut; reuse accepted ppRef fits."""
from pathlib import Path
import json,math,hashlib,shutil
import numpy as np
import uproot
import ROOT as R
repo=Path(__file__).resolve().parents[2]
out=repo/'fitER/results/ppref/common_precut_mc'
assert not out.exists(), 'Refusing to overwrite results'
ref=repo/'fitER/results/ppref/two_peak_precut'
rref=json.loads((ref/'fit_result.json').read_text())
assert json.loads((ref/'audit.json').read_text())['status']=='PASS'
selection='Bpt > 7.5 && Bpt < 50 && abs(By) < 2.4 && BQvalue < 0.15'
assert json.loads((ref/'manifest.json').read_text())['selection']==selection
out.mkdir(parents=True);art=out/'artifacts';art.mkdir();shutil.copy2(__file__,art)
shutil.copy2(repo/'fitER/models/TwoPeakPlots.C',art/'TwoPeakPlots.C')
R.gROOT.SetBatch(True);R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
R.gSystem.SetBuildDir(str(art),True);assert R.gSystem.CompileMacro(str(art/'TwoPeakPlots.C'),'k')>0
result={'selection':selection,'ml_cut':None,'weight':'pThatreweight','sample':'prompt RECO PbPb23','mc':{},'inputs':{},'fits':{},'reference':{'path':str(ref),'sha256':hashlib.sha256((ref/'fit_result.json').read_bytes()).hexdigest()},'configuration':{'mean_half_range':.01,'sigma_initial':[.01,.005],'sigma_range':[.001,.1],'fraction_initial':.5,'fraction_range':[.01,1.],'scale':1,'strategy':2,'hesse':True,'SumW2Error':False}}
keep=[]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
def record(q):return {'status':q.status(),'covQual':q.covQual(),'edm':q.edm(),'nll':q.minNll(),'parameters':{v.GetName():{'value':v.getVal(),'error':v.getError(),'range':[v.getMin(),v.getMax()]} for v in q.floatParsFinal()}}
def passed(q):return q.status()==0 and q.covQual()==3 and math.isfinite(q.edm()) and q.edm()<.001
for key,particle,nom,lo,hi in [('x','X3872',3.87169,3.84,3.9),('psi','PSI2S',3.68610,3.65610,3.71610)]:
 src=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23')/('flat_ntmix_PbPb23_MC_'+particle+'.root');tree='ntmix_'+particle;st=src.stat()
 with uproot.open(src) as f:a=f[tree].arrays(['Bmass','Bpt','By','BQvalue','pThatreweight'],library='np')
 fid=(a['Bpt']>7.5)&(a['Bpt']<50)&(abs(a['By'])<2.4)&(a['BQvalue']<.15)
 mask=fid&(a['Bmass']>3.62)&(a['Bmass']<4.0)
 b={n:a[n][mask] for n in ['Bmass','pThatreweight']};b['source_entry']=np.flatnonzero(mask).astype('int64');weights=b['pThatreweight'];assert np.isfinite(weights).all() and (weights>0).all()
 cache=art/(key+'_cache.root')
 with uproot.recreate(cache) as f:f.mktree(tree,{n:v.dtype for n,v in b.items()});f[tree].extend(b)
 with uproot.open(cache) as f:check=f[tree].arrays(library='np')
 for n in b:assert np.array_equal(b[n],check[n])
 mass=R.RooRealVar('Bmass','Bmass',3.62,4.0);mass.setRange('mc_range',lo,hi);weight=R.RooRealVar('pThatreweight','pThatreweight',-1e6,1e6)
 ff=R.TFile.Open(str(cache));d=R.RooDataSet(key+'_mc','',R.RooArgSet(mass,weight),R.RooFit.Import(ff.Get(tree)),R.RooFit.WeightVar('pThatreweight'))
 assert d.numEntries()==len(weights);assert abs(d.sumEntries()-weights.sum())<1e-6*max(1,weights.sum())
 mu=R.RooRealVar(key+'_mc_mean','mean',nom,nom-.01,nom+.01);s1=R.RooRealVar(key+'_sigma1','sigma1',.01,.001,.1);s2=R.RooRealVar(key+'_sigma2','sigma2',.005,.001,.1);fraction=R.RooRealVar(key+'_fraction','fraction',.5,.01,1.)
 g1=R.RooGaussian(key+'_g1','',mass,mu,s1);g2=R.RooGaussian(key+'_g2','',mass,mu,s2);pdf=R.RooAddPdf(key+'_mc_pdf','',R.RooArgList(g1,g2),R.RooArgList(fraction))
 opts=[R.RooFit.Save(),R.RooFit.Range('mc_range'),R.RooFit.Strategy(2),R.RooFit.Hesse(True),R.RooFit.PrintLevel(-1),R.RooFit.Warnings(False),R.RooFit.Verbose(False),R.RooFit.SumW2Error(False)]
 fit=pdf.fitTo(d,*opts);attempts=[record(fit)];first=None
 if math.isfinite(fit.minNll()) and not passed(fit):
  first=fit;pdf.getParameters(d).assignValueOnly(fit.floatParsFinal());fit=pdf.fitTo(d,*opts);attempts.append(record(fit))
 nll=pdf.createNLL(d,R.RooFit.Range('mc_range'));assert abs(nll.getVal()-fit.minNll())<1e-6
 values={'mean':mu.getVal(),'sigma1':s1.getVal(),'sigma2':s2.getVal(),'fraction':fraction.getVal(),'range':[lo,hi]}
 values['sigma_eff']=math.sqrt(fraction.getVal()*s1.getVal()**2+(1-fraction.getVal())*s2.getVal()**2)
 values['parameter_boundary']=any(min(v.getVal()-v.getMin(),v.getMax()-v.getVal())<=1e-4*(v.getMax()-v.getMin()) for v in [mu,s1,s2,fraction])
 result['mc'][key]=values;result['fits'][key]={'attempts':attempts,'quality_pass':passed(fit),'nll_recomputed':True}
 win=(b['Bmass']>lo)&(b['Bmass']<hi);ww=weights[win]
 result['inputs'][key]={'path':str(src),'tree':tree,'size':st.st_size,'mtime_ns':st.st_mtime_ns,'sha256':sha(src),'source_entries':len(mask),'precut_entries':int(fid.sum()),'cache_entries':len(weights),'fit_window_entries':int(win.sum()),'fit_window_sumw':float(ww.sum()),'fit_window_neff':float(ww.sum()**2/(ww@ww)),'weight_min':float(weights.min()),'weight_max':float(weights.max()),'cache':str(cache),'cache_sha256':sha(cache),'cache_roundtrip_exact':True}
 assert src.stat().st_mtime_ns==st.st_mtime_ns
 output=R.TFile.Open(str(out/(key+'_fit_workspace.root')),'RECREATE');w=R.RooWorkspace('ws_mc');imp=getattr(w,'import');imp(d);imp(pdf);w.Write();fit.Write('fit_result_mc')
 if first:first.Write('fit_result_mc_first_attempt')
 output.Close()
 R.PlotPeakMC(str(out/(key+'_mc_template_fit.pdf')),d,pdf,mass,fit,mu.getVal(),s1.getVal(),s2.getVal(),fraction.getVal(),lo,hi,'PbPb X' if key=='x' else 'PbPb #psi(2S)')
 keep.extend([ff,d,mass,weight,mu,s1,s2,fraction,g1,g2,pdf,fit,first,nll])
 (out/'fit_result.json').write_text(json.dumps(result,indent=2)+'\n')
result['quality_pass']=all(x['quality_pass'] for x in result['fits'].values());(out/'fit_result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2),flush=True)
assert result['quality_pass'],'MC fit quality failed; review before comparison'
