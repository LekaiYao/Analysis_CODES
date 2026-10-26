#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""Authorized v27/v1 xeff50 three-component sPlot, with frozen accepted shapes."""
from pathlib import Path
import argparse,json,shutil,math,re,hashlib
import numpy as np
import uproot,numexpr,ROOT as R
from fit_input_contract import load_training,selection_expression
from ppref_psi2s_pthat_validation import fingerprint,dump,scalar,stats,record

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=args.output.resolve();assert not out.exists()
 repo=Path(__file__).resolve().parents[2];tag='X_pb23_v27_fid18_9v9_rw0_xgb_v1';base=repo/'fitER/results/ml_fits'/tag;source=base/'checks/two_peak_mc_width/xeff50'
 R.gROOT.SetBatch(True);assert R.gROOT.GetVersion()=='6.32.02';R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
 training,points=load_training(repo,tag);ref=json.loads((source/'fit_result.json').read_text());manifest=json.loads((source/'manifest.json').read_text());parent=base/'two_peak/xeff50/manifest.json';pm=json.loads(parent.read_text());cut=manifest['threshold']
 assert ref['audit']['quality_pass'] and training['selection']==manifest['selection']==pm['selection']
 assert next(p['score_threshold'] for p in points if abs(p['target_efficiency']-.5)<1e-9)==cut==pm['threshold']
 assert fingerprint(parent)['sha256']==manifest['source_manifest_sha256'];assert fingerprint(Path(manifest['source_workspace']))['sha256']==manifest['source_sha256']
 paths={k:Path(v['source']) for k,v in pm['inputs'].items()};protected={**paths,'workspace':source/'fit_workspace.root','fit_result':source/'fit_result.json','manifest':source/'manifest.json','parent_manifest':parent,'thresholds':Path(training['thresholds_file']),'plot_config':repo/'plotER/validation_archive/aux.h','script':Path(__file__)}
 provenance={k:fingerprint(p) for k,p in protected.items()};names={}
 for k,p in paths.items():
  st=p.stat();old=pm['inputs'][k];assert (st.st_size,st.st_mtime_ns)==(old['size'],old['mtime_ns']),k
  with uproot.open(p) as f:names[k]={n for n in f[old['tree']].keys() if scalar(f[old['tree']][n])}
 variables=sorted(set.intersection(*names.values())-{'Bmass','Bgen','pThatreweight'});assert len(variables)==68 and 'Prediction' in variables
 arrays={};counts={}
 for k,p in paths.items():
  columns=sorted(set(variables)|{'Bmass'}|({'pThatreweight'} if k!='data' else set()))
  with uproot.open(p) as f:a=f[pm['inputs'][k]['tree']].arrays(columns,library='np')
  fid=numexpr.evaluate(selection_expression(training['selection']),local_dict=a);assert fid.all();assert np.isfinite(a['Bmass']).all() and np.isfinite(a['Prediction']).all()
  mask=fid&(a['Prediction']>cut)
  if k=='data':mask &= (a['Bmass']>3.62)&(a['Bmass']<4.0)
  arrays[k]={n:a[n][mask] for n in columns};arrays[k]['source_entry']=np.flatnonzero(mask).astype('int64');counts[k]={'source':len(mask),'selected':int(mask.sum())}
  old=pm['inputs'][k];assert fingerprint(Path(old['cache']))['sha256']==old['sha256']
  with uproot.open(old['cache']) as f:b=f[old['tree']].arrays(library='np')
  full=arrays[k];massmask=(full['Bmass']>3.62)&(full['Bmass']<4.0)
  for name in b:assert np.array_equal(full[name][massmask],b[name]),(k,name)
  if k!='data':assert np.isfinite(full['pThatreweight']).all() and np.all(full['pThatreweight']>0)
 d=arrays['data'];assert len(d['Bmass'])==7017
 out.mkdir(parents=True);art=out/'artifacts';art.mkdir();dump(out/'provenance.json',provenance)
 for p in [Path(__file__),repo/'scripts/workflows/fit_input_contract.py',repo/'scripts/workflows/ppref_psi2s_pthat_validation.py',repo/'plotER/validation_archive/aux.h']:shutil.copy2(p,art/p.name)
 inventory=[{'variable':n,'included':n in variables,'reason':'common_RECO_scalar_or_score' if n in variables else 'mass_discriminant_or_MC_only_weight_truth'} for n in sorted(set.union(*names.values()))]
 inventory += [{'variable':n,'included':False,'reason':'missing_no_alias'} for n in ['Bmu1y','Bmu2y'] if n not in set.union(*names.values())];dump(out/'variable_inventory.json',inventory)
 f=R.TFile.Open(str(source/'fit_workspace.root'));w=f.Get('ws_nominal');model=w.pdf('model');data=w.data('data_data');mass=w.var('Bmass');orig=f.Get('fit_result_alt');pars=model.getParameters(data);pars.assignValueOnly(orig.floatParsFinal());assert data.numEntries()==len(d['Bmass'])
 for i,value in enumerate(d['Bmass']):assert data.get(i).getRealValue('Bmass')==value
 ykeys=['x_yield','psi_yield','nbkg'];yields=[w.var(k) for k in ykeys]
 for p in pars:p.setConstant(p.GetName() not in ykeys)
 fixed={p.GetName():p.getVal() for p in pars if p.isConstant()};assert fixed['x_scale']==fixed['psi_scale']==1.
 for name,val in manifest['fixed_parameters'].items():assert fixed[name]==val
 opts=[R.RooFit.Save(),R.RooFit.Extended(True),R.RooFit.Range('all'),R.RooFit.Strategy(2),R.RooFit.Hesse(True),R.RooFit.PrintLevel(-1)]
 good=lambda q:q.status()==0 and q.covQual()==3 and math.isfinite(q.edm()) and q.edm()<1e-3
 fit=model.fitTo(data,*opts);history=[record(fit)]
 if not good(fit) and math.isfinite(fit.minNll()):pars.assignValueOnly(fit.floatParsFinal());fit=model.fitTo(data,*opts);history.append(record(fit))
 dump(out/'retry_history.json',history);assert good(fit),history;assert set(p.GetName() for p in fit.floatParsFinal())==set(ykeys)
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
 summary={'status':'SPLOT_PASS','tag':tag,'xeff':50,'score_threshold':cut,'selection':training['selection'],'data_mass_range':[3.62,4.0],'mc_validation_mass_cut':None,'mc_weight':'pThatreweight','fit_source':str(source),'fit':record(fit),'splot_yields':dict(zip(ykeys,nv.tolist())),'sweights':results,'mc':{k:stats(arrays[k]['pThatreweight'].astype(float)) for k in ['x','psi']},'cutflow':counts,'variables':variables,'fixed_parameters':fixed,'checks':{'source_entry_alignment':'PASS','analytic_sweight_max_difference':float(np.max(abs(sw-analytic))),'sumw_minus_yield':(sw.sum(axis=0)-nv).tolist(),'event_sum_max_difference':float(np.max(abs(sw.sum(axis=1)-1))),'information_condition':float(np.linalg.cond(information))},'limitations':['Inclusive DATA signal versus prompt MC, not prompt-only closure.','Signed sWeights retained; mass independence assumed after score selection.','Low effective signal statistics; CDF distances are descriptive, not KS p-values.','Fixed shape uncertainty not propagated; no bootstrap or retraining.']}
 dump(out/'summary.json',summary)
 for k,p in protected.items():assert fingerprint(p)==provenance[k],k
 f.Close();print('SPLOT COMPLETE',json.dumps(results),flush=True)
if __name__=='__main__':main()
