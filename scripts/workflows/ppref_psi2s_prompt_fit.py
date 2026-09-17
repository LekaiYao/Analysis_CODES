#!/usr/bin/env python3
"""Prompt-enriched ppRef Psi2S fit/sPlot: baseline selection plus SV-PV significance < 2.
Run in ROOT 6.32.02 LCG_106. Existing output directory is never overwritten.
"""
import argparse, csv, hashlib, json, math
from pathlib import Path
import numpy as np
import uproot
import ROOT as R
from scipy.special import ndtr

BASE=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24')
FILES={'data':BASE/'flat_ntmix_ppRef_DATA.root','mc':BASE/'flat_ntmix_ppRef_MC_PSI2S.root'}
TREES={'data':'ntmix','mc':'ntmix_PSI2S'}
EXPECTED={'data':'03706f676cf24d3bf0c24c5ab49cf54345369c8cfb80d6f75d3d287e5787e599','mc':'8777f8c64135aa5e0b829c6e21796aa7619be16d8f4d3ddec813553eb3c8d837'}
SEL='Bpt > 7.5 && Bpt < 50 && abs(By) < 2.4 && BQvalue < 0.15 && Bnorm_svpvDistance_2D < 2'
M0=3.68610

def fingerprint(p):
 a=p.stat();h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 z=p.stat();assert (a.st_size,a.st_mtime_ns)==(z.st_size,z.st_mtime_ns)
 return {'path':str(p),'bytes':z.st_size,'mtime_ns':z.st_mtime_ns,'sha256':h.hexdigest()}

def dump(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def scalar(b):
 try:return np.dtype(b.interpretation.numpy_dtype).shape==() and np.dtype(b.interpretation.numpy_dtype).kind in 'biuf'
 except Exception:return False

def select(a):return (a['Bpt']>7.5)&(a['Bpt']<50)&(np.abs(a['By'])<2.4)&(a['BQvalue']<.15)&(a['Bnorm_svpvDistance_2D']<2)
def stats(w):return {'entries':len(w),'sumw':float(w.sum()),'sumw2':float(w@w),'neff':float(w.sum()**2/(w@w)),'negative_fraction':float(np.mean(w<0))}
def record(f):
 pars={}
 for p in f.floatParsFinal():
  pars[p.GetName()]={'value':p.getVal(),'error':p.getError(),'min':p.getMin(),'max':p.getMax(),'near_boundary':min(p.getVal()-p.getMin(),p.getMax()-p.getVal())<1e-4*(p.getMax()-p.getMin())}
 return {'status':f.status(),'covQual':f.covQual(),'edm':f.edm(),'minNll':f.minNll(),'parameters':pars}

def fit(pdf,ds,weighted=False):
 opts=[R.RooFit.Save(),R.RooFit.PrintLevel(-1),R.RooFit.Strategy(2),R.RooFit.Offset(True)]
 if weighted:opts.append(R.RooFit.SumW2Error(True))
 f=pdf.fitTo(ds,*opts)
 if f.status()!=0 or f.covQual()!=3 or f.edm()>1e-3:f=pdf.fitTo(ds,*opts)
 if f.status()!=0 or f.covQual()!=3 or f.edm()>1e-3:raise RuntimeError('fit quality failed: '+str(record(f)))
 return f

def cdf(x,w,y,v):
 ix=np.argsort(x,kind='stable');iy=np.argsort(y,kind='stable');z=np.union1d(x,y)
 cx=np.r_[0.,np.cumsum(w[ix])][np.searchsorted(x[ix],z,side='right')]/w.sum()
 cy=np.r_[0.,np.cumsum(v[iy])][np.searchsorted(y[iy],z,side='right')]/v.sum()
 return z,cx,cy,float(np.max(np.abs(cx-cy)))

def corr(x,y,w):
 sx=x-np.average(x,weights=w);sy=y-np.average(y,weights=w);den=np.sqrt(np.sum(w*sx*sx)*np.sum(w*sy*sy))
 return float(np.sum(w*sx*sy)/den) if den>0 else None

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args();out=Path(args.output)
 out.mkdir(parents=True,exist_ok=False)
 for d in ['figures','artifacts']: (out/d).mkdir()
 R.gROOT.SetBatch(True)
 assert R.gROOT.GetVersion()=='6.32.02'
 R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
 provenance={k:fingerprint(p) for k,p in FILES.items()}
 for k in FILES:assert provenance[k]['sha256']==EXPECTED[k],f'{k} input changed from inspected version'
 provenance['script']=fingerprint(Path(__file__))
 dump(out/'provenance.json',provenance)
 inventory=[];names={}
 for k,p in FILES.items():
  with uproot.open(p) as f:names[k]={n for n in f[TREES[k]].keys() if scalar(f[TREES[k]][n])}
 variables=sorted((names['data']&names['mc'])-{'Bmass','pThatreweight','Bgen'})
 for n in sorted(names['data']|names['mc']):inventory.append({'variable':n,'included':n in variables,'reason':'common_RECO_scalar' if n in variables else 'mass_discriminant_or_noncommon_or_weight'})
 for n in ['Bmu1y','Bmu2y']:
  if n not in names['data']|names['mc']:inventory.append({'variable':n,'included':False,'reason':'missing_no_alias'})
 dump(out/'variable_inventory.json',inventory)
 arrays={};cutflow={}
 for k,p in FILES.items():
  read=sorted(set(variables)|{'Bmass','Bpt','By','BQvalue','Bnorm_svpvDistance_2D'}|({'pThatreweight'} if k=='mc' else set()))
  chunks={n:[] for n in read};chunks['source_entry']=[];offset=0;nbase=0;total=0
  with uproot.open(p) as f:
   for a in f[TREES[k]].iterate(read,step_size=100000,library='np'):
    mask=select(a);nbase+=int(mask.sum());n=len(mask);total+=n
    if k=='data':mask &= (a['Bmass']>3.6)&(a['Bmass']<3.8)
    for key in read:chunks[key].append(a[key][mask])
    chunks['source_entry'].append(np.arange(offset,offset+n,dtype=np.int64)[mask]);offset+=n
  arrays[k]={n:np.concatenate(v) for n,v in chunks.items()};cutflow[k]={'total_RECO':total,'preselection':nbase,'retained':len(arrays[k]['Bmass'])}
 d=arrays['data'];mc=arrays['mc'];w=mc['pThatreweight'].astype(float)
 assert np.isfinite(w).all() and (w>0).all()
 assert np.isfinite(d['Bmass']).all() and np.isfinite(mc['Bmass']).all()
 print('INPUT',cutflow,flush=True)
 mass=R.RooRealVar('Bmass','Bmass',3.6,3.8);obs=R.RooArgSet(mass)
 peak=(mc['Bmass']>M0-.03)&(mc['Bmass']<M0+.03)
 mass.setRange('mc_peak',M0-.03,M0+.03)
 wt=R.RooRealVar('pThatreweight','pThatreweight',0.,float(w.max()*2))
 md=R.RooDataSet.from_numpy({'Bmass':mc['Bmass'][peak],'pThatreweight':w[peak]},[mass,wt],name='weighted_mc',weight_name='pThatreweight')
 data=R.RooDataSet.from_numpy({'Bmass':d['Bmass']},[mass],name='data')
 mean=R.RooRealVar('mean','mean',M0,M0-.01,M0+.01)
 s1=R.RooRealVar('sigma1','sigma1',.01,.001,.1);s2=R.RooRealVar('sigma2','sigma2',.005,.001,.1)
 frac=R.RooRealVar('fraction','fraction',.5,.01,1.)
 scale=R.RooRealVar('scale','scale',1.,.9,1.15);scale.setConstant(True)
 a1=R.RooFormulaVar('scaled_sigma1','@0*@1',R.RooArgList(s1,scale));a2=R.RooFormulaVar('scaled_sigma2','@0*@1',R.RooArgList(s2,scale))
 g1=R.RooGaussian('g1','g1',mass,mean,a1);g2=R.RooGaussian('g2','g2',mass,mean,a2)
 sig=R.RooAddPdf('signal','signal',R.RooArgList(g1,g2),R.RooArgList(frac))
 # Same nominal MC peak interval; weighted shape likelihood with SumW2 covariance.
 mcopts=[R.RooFit.Range('mc_peak'),R.RooFit.Save(),R.RooFit.Strategy(2),R.RooFit.PrintLevel(-1),R.RooFit.Offset(True)]
 # SumW2Error's auxiliary squared-weight Hesse can report an EDM unrelated
 # to minimization of the original weighted likelihood. Audit both explicitly.
 fmc_raw=sig.fitTo(md,*mcopts,R.RooFit.SumW2Error(False))
 if fmc_raw.status()!=0 or fmc_raw.covQual()!=3 or fmc_raw.edm()>1e-3:
  fmc_raw=sig.fitTo(md,*mcopts,R.RooFit.SumW2Error(False))
 dump(out/'mc_fit_raw_likelihood.json',record(fmc_raw))
 assert fmc_raw.status()==0 and fmc_raw.covQual()==3 and fmc_raw.edm()<1e-3,record(fmc_raw)
 raw_mc_values={p.GetName():p.getVal() for p in fmc_raw.floatParsFinal()}
 fmc=sig.fitTo(md,*mcopts,R.RooFit.SumW2Error(True))
 mcr=record(fmc);mcr['original_likelihood_edm']=fmc_raw.edm()
 mcr['sumw2_parameter_shift']={p.GetName():p.getVal()-raw_mc_values[p.GetName()] for p in fmc.floatParsFinal()}
 dump(out/'mc_fit.json',mcr)
 assert fmc.status()==0 and fmc.covQual()==3,mcr
 for p in fmc.floatParsFinal():assert abs(mcr['sumw2_parameter_shift'][p.GetName()])<.1*p.getError(),mcr
 mcpars=(mean.getVal(),s1.getVal(),s2.getVal(),frac.getVal())
 for x in [s1,s2,frac]:x.setConstant(True)
 scale.setConstant(False)
 c0=R.RooRealVar('c0','c0',-.35,-2.,2.);c1=R.RooRealVar('c1','c1',-.05,-2.,2.)
 bkg=R.RooChebychev('background','background',mass,R.RooArgList(c0,c1))
 near=int(np.sum(np.abs(d['Bmass']-M0)<.005));n=len(d['Bmass'])
 ns=R.RooRealVar('nsig','nsig',near*.4,0.,max(10.,near*2.));nb=R.RooRealVar('nbkg','nbkg',n*.7,n*.1,float(n))
 model=R.RooAddPdf('model','model',R.RooArgList(sig,bkg),R.RooArgList(ns,nb))
 fd=fit(model,data);fr=record(fd);dump(out/'data_fit.json',fr)
 # Check exact quadratic minimum, including interior extremum, before sPlot.
 z=[-1.,1.]
 if c1.getVal()!=0:
  vertex=-c0.getVal()/(4*c1.getVal())
  if -1<vertex<1:z.append(vertex)
 assert min(1+c0.getVal()*q+c1.getVal()*(2*q*q-1) for q in z)>0
 for x in [mean,scale,c0,c1]:x.setConstant(True)
 fs=fit(model,data);sr=record(fs)
 swobj=R.RooStats.SPlot('splot','splot',data,model,R.RooArgList(ns,nb))
 sw=np.fromiter((swobj.GetSWeight(i,'nsig') for i in range(n)),dtype=float,count=n)
 bw=np.fromiter((swobj.GetSWeight(i,'nbkg') for i in range(n)),dtype=float,count=n)
 assert np.isfinite(sw).all() and np.isfinite(bw).all()
 assert abs(sw.sum()-ns.getVal())<1e-5*max(1.,ns.getVal())
 assert abs(bw.sum()-nb.getVal())<1e-5*max(1.,nb.getVal())
 assert np.max(np.abs(sw+bw-1))<1e-3
 d['signal_sWeight']=sw
 with uproot.recreate(out/'artifacts/psi2s_sweights.root') as f:
  f.mktree('ntmix_PSI2S_sWeight',{k:v.dtype for k,v in d.items()});f['ntmix_PSI2S_sWeight'].extend(d)
 with uproot.recreate(out/'artifacts/selected_mc.root') as f:
  f.mktree('ntmix_PSI2S',{k:v.dtype for k,v in mc.items()});f['ntmix_PSI2S'].extend(mc)
 rf=R.TFile(str(out/'artifacts/fit_workspace.root'),'RECREATE');ws=R.RooWorkspace('ws_nominal')
 getattr(ws,'import')(model);getattr(ws,'import')(data);ws.Write();fmc.Write('fit_mc');fd.Write('fit_data');fs.Write('fit_splot');rf.Close()
 summary={'selection':SEL,'data_mass_range':[3.6,3.8],'mc_validation_mass_cut':None,'mc_shape_range':[M0-.03,M0+.03],'mc_weight':'pThatreweight','cutflow':cutflow,'mc':stats(w),'splot':stats(sw),'data_fit':fr,'mc_fit':mcr,'splot_fit':sr,'sumw_minus_yield':float(sw.sum()-ns.getVal()),'max_event_sweight_sum_deviation':float(np.max(np.abs(sw+bw-1))),'common_variables':len(variables),'warnings':['Prompt-enriched DATA compared with prompt MC; residual nonprompt is not subtracted.','sPlot assumes mass independence; CDF distances are descriptive, not KS p-values.','No bootstrap, ML training, automatic variable choice, or significance calculation.']}
 if any(v['near_boundary'] for v in fr['parameters'].values()):summary['warnings'].append('DATA fit parameter near nominal boundary; no range tuning performed.')
 dump(out/'summary.json',summary)
 for k,p in FILES.items():assert fingerprint(p)==provenance[k]
 assert fingerprint(Path(__file__))==provenance['script']
 summary['status']='fit_and_splot_complete';dump(out/'summary.json',summary)
 print('FIT/SPLOT COMPLETE',out,summary['splot'],flush=True)
if __name__=='__main__':main()
