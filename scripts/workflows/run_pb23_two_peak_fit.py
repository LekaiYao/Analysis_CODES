from pathlib import Path
import json,hashlib,shutil,math,re,argparse,subprocess,sys
from fit_input_contract import load_training, selection_expression
ap=argparse.ArgumentParser(description='Explicitly requested single-point X+psi(2S) fit')
ap.add_argument('--thresholds',type=Path,help='Validated downstream threshold supplement for an explicitly selected point');ap.add_argument('--tag');ap.add_argument('--xeff',type=int);ap.add_argument('--ppref',action='store_true',help='ppRef ordinary preselection, no ML cut, separate results/report');ap.add_argument('--float-width',action='store_true',help='Explicit legacy PbPb strategy; ppRef keeps its existing floating-width model');ap.add_argument('--dry-run',action='store_true');args=ap.parse_args()
repo=Path(__file__).resolve().parents[2];fixed_width=not args.ppref and not args.float_width
if args.ppref:
 assert args.tag is None and args.xeff is None, 'ppRef has no ML tag or xeff cut'
 tag='ppRef_precut';out=repo/'fitER/results/ppref/two_peak_precut'
else:
 tag=args.tag;assert tag and re.fullmatch(r'X_pb23_[A-Za-z0-9_]+',tag);assert args.xeff is not None and 0<args.xeff<100
 out=repo/'fitER/results/ml_fits'/tag/('two_peak_mc_width' if fixed_width else 'two_peak')/f'xeff{args.xeff}'
if args.dry_run:print(out,'mass=[3.62,4.00]; independent means; '+('both scales fixed to 1' if fixed_width else 'independent floating scales')+'; then conditional test');raise SystemExit(0)
assert not out.exists(), 'Refusing to overwrite existing result'
if args.ppref:
 base=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24');thr=None
 t={'selection':'Bpt > 7.5 && Bpt < 50 && abs(By) < 2.4 && BQvalue < 0.15','weight':'pThatreweight','thresholds_file':None,'thresholds_sha256':None}
else:
 t,points=load_training(repo,tag)
 if args.thresholds:
  extra=json.loads(args.thresholds.read_text());assert extra['schema_version']==1 and extra['train_tag']==tag
  assert extra['parent_thresholds_sha256']==t['thresholds_sha256'] and extra['selection']==t['selection']
  assert extra['weight_branch']==[t['weight']] and extra['comparison']=='Prediction > threshold'
  assert extra['existing_thresholds_exactly_reproduced']==len(points)
  assert hashlib.sha256(Path(extra['reference_file']).read_bytes()).hexdigest()==extra['reference_sha256']
  for p in points:
   for q in extra['thresholds']:
    if abs(p['target_efficiency']-q['target_efficiency'])<1e-8:assert p['score_threshold']==q['score_threshold']
  t['thresholds_file']=str(args.thresholds.resolve());t['thresholds_sha256']=hashlib.sha256(args.thresholds.read_bytes()).hexdigest();points=extra['thresholds']
 matches=[p for p in points if abs(p['target_efficiency']-args.xeff/100)<1e-8];assert len(matches)==1;thr=matches[0]['score_threshold']
 assert t['weight']=='pThatreweight', 'This validated two-peak contract requires pThatreweight for both MC samples'
 base=repo.parent/'XGBoost/output/selected'/tag
 assert (base/'MC_psi2s_with_score.root').is_file(), 'Matching scored psi(2S) MC is required'
import numpy as np,uproot,numexpr,ROOT as R
R.gROOT.SetBatch(True);R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
out.mkdir(parents=True);art=out/'artifacts';art.mkdir();shutil.copy2(__file__,art);shutil.copy2(repo/'scripts/workflows/fit_input_contract.py',art)
records={};arrays={};handles=[];keep=[]
for key,name,tree in [('data','DATA_with_score.root','ntmix'),('x','MC_with_score.root','ntmix_X3872'),('psi','MC_psi2s_with_score.root','ntmix_PSI2S')]:
 if args.ppref:
  src=base/('flat_ntmix_ppRef_'+{'data':'DATA','x':'MC_X3872','psi':'MC_PSI2S'}[key]+'.root')
 else:
  src=Path(t['sources']['data' if key=='data' else 'mc']['path']) if key!='psi' else base/name;tree=t['sources']['data' if key=='data' else 'mc']['tree'] if key!='psi' else tree
 stat=src.stat();kept=['Bmass']+([] if args.ppref else ['Prediction'])+(['pThatreweight'] if key!='data' else [])
 branches=sorted(set(kept+re.findall(r'\bB\w+\b',t['selection'])))
 with uproot.open(src) as f:a=f[tree].arrays(branches,library='np')
 fid=numexpr.evaluate(selection_expression(t['selection']),local_dict=a)
 if not args.ppref:assert fid.all(),key
 mask=fid&(True if args.ppref else a['Prediction']>thr)&(a['Bmass']>3.62)&(a['Bmass']<4.0);a={n:a[n][mask] for n in kept};a['source_entry']=np.flatnonzero(mask).astype('int64');cache=art/(key+'_cache.root')
 with uproot.recreate(cache) as f:f.mktree(tree,{n:v.dtype for n,v in a.items()});f[tree].extend(a)
 assert src.stat().st_mtime_ns==stat.st_mtime_ns
 records[key]={'source':str(src),'size':stat.st_size,'mtime_ns':stat.st_mtime_ns,'source_entries':len(mask),'entries':int(mask.sum()),'cache':str(cache),'sha256':hashlib.sha256(cache.read_bytes()).hexdigest(),'tree':tree};arrays[key]=a
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
 mcresults[key]={'mean':mu.getVal(),'sigma1':s1.getVal(),'sigma2':s2.getVal(),'fraction':frac.getVal(),'range':[lo,hi],'sumw':float(arrays[key]['pThatreweight'].sum()),'sumw2':float((arrays[key]['pThatreweight']**2).sum())}
 for p in [mu,s1,s2,frac]:p.setConstant(True)
 mean=var(key+'_mean',mu.getVal(),nom-.01,nom+.01);scale=var(key+'_scale',1.,.5,1.5);ss1=R.RooProduct(key+'_scaled_sigma1','',R.RooArgList(scale,s1));ss2=R.RooProduct(key+'_scaled_sigma2','',R.RooArgList(scale,s2));dg1=R.RooGaussian(key+'_g1','',mass,mean,ss1);dg2=R.RooGaussian(key+'_g2','',mass,mean,ss2);sig=R.RooAddPdf(key+'_signal','',R.RooArgList(dg1,dg2),R.RooArgList(frac));near=int(np.sum(abs(arrays['data']['Bmass']-nom)<.005));yieldv=var(key+'_yield',max(0.,.4*near),0.,max(10.,2.*near));params[key]={'mean':mean,'scale':scale,'yield':yieldv};pdfs[key]=sig;keep.extend([ss1,ss2,dg1,dg2,sig]);pdfs[key+'_mc']=pdf
if fixed_width:
 for pp in params.values():pp['scale'].setVal(1.);pp['scale'].setConstant(True)
n=data.numEntries();a0=var('a0',-.35,-2,2);a1=var('a1',-.05,-2,2);nb=var('nbkg',.7*n,.1*n,n);bg=R.RooChebychev('backgroundPdf','',mass,R.RooArgList(a0,a1));model=R.RooAddPdf('model','',R.RooArgList(pdfs['x'],pdfs['psi'],bg),R.RooArgList(params['x']['yield'],params['psi']['yield'],nb));keep.extend([bg,model]);alt=fit(model,data,'all');fitresults['alt']=alt;allpars=model.getParameters(data);snapshot=allpars.snapshot()
def quality(r):return {'status':r.status(),'covQual':r.covQual(),'edm':r.edm(),'nll':r.minNll()}
def boundary(p):return min(abs(p.getVal()-p.getMin()),abs(p.getVal()-p.getMax()))<=1e-4*(p.getMax()-p.getMin())
result={'tag':tag,'xeff':args.xeff,'threshold':thr,'mass_range':[3.62,4.0],'mc':mcresults,'data_entries':n,'fits':{'alt':quality(alt)},'signals':{k:{n:{'value':p.getVal(),'error':p.getError(),'range':[p.getMin(),p.getMax()],'at_boundary':boundary(p)} for n,p in pp.items()} for k,pp in params.items()},'background':{p.GetName():{'value':p.getVal(),'error':p.getError(),'at_boundary':boundary(p)} for p in [a0,a1,nb]}}
for key in ['x','psi']:
 allpars.assignValueOnly(snapshot)
 for pp in params.values():
  for name,p in pp.items():p.setConstant(fixed_width and name=='scale')
 params[key]['yield'].setVal(0.)
 for p in params[key].values():p.setConstant(True)
 rr=fit(model,data,'all');fitresults['null_'+key]=rr;result['fits']['null_'+key]=quality(rr);q=max(0.,2*(rr.minNll()-alt.minNll()));result['signals'][key].update(q0=q,Z_PL=math.sqrt(q))
 other='psi' if key=='x' else 'x';assert {p.GetName() for p in rr.floatParsFinal()}==({other+'_mean',other+'_yield','a0','a1','nbkg'} | (set() if fixed_width else {other+'_scale'}))
 check_nll=model.createNLL(data,R.RooFit.Extended(True),R.RooFit.Range('all'));assert abs(check_nll.getVal()-rr.minNll())<1e-5
allpars.assignValueOnly(snapshot)
for pp in params.values():
 for name,p in pp.items():p.setConstant(fixed_width and name=='scale')
for k in ['x','psi']:result['fits'][k+'_mc']=quality(fitresults[k+'_mc'])
assert {p.GetName() for p in alt.floatParsFinal()}==({'x_mean','x_yield','psi_mean','psi_yield','a0','a1','nbkg'} | (set() if fixed_width else {'x_scale','psi_scale'}))
if fixed_width:
 for k in ['x','psi']:
  assert params[k]['scale'].isConstant() and params[k]['scale'].getVal()==1.
  result['signals'][k]['scale'].update(fixed=True,error=None)
result['shape_strategy']='mc_width_mean_float' if fixed_width else 'mean_scale_float'
for k in ['x','psi']:assert abs(mc[k].sumEntries()-mcresults[k]['sumw'])<1e-5*max(1,mcresults[k]['sumw'])
pts=[-1.,1.];aa=a0.getVal();bb=a1.getVal()
if bb and -1 < -aa/(4*bb) < 1:pts.append(-aa/(4*bb))
positive=min(1+aa*x+bb*(2*x*x-1) for x in pts);assert positive>0
nll=model.createNLL(data,R.RooFit.Extended(True),R.RooFit.Range('all'));assert abs(nll.getVal()-alt.minNll())<1e-5
result['audit']={'quality_pass':all(q['status']==0 and q['covQual']==3 and q['edm']<1e-3 for q in result['fits'].values()),'counts_weights_parameters_nll':'PASS','background_min':positive}
for info in records.values():assert hashlib.sha256(Path(info['cache']).read_bytes()).hexdigest()==info['sha256']
f=R.TFile.Open(str(out/'fit_workspace.root'),'RECREATE');w=R.RooWorkspace('ws_nominal');imp=getattr(w,'import');imp(data);imp(model)
for k in ['x','psi']:imp(mc[k]);imp(pdfs[k+'_mc'],R.RooFit.RecycleConflictNodes())
w.Write()
for k,r in fitresults.items():r.Write('fit_result_'+k)
f.Close();(out/'fit_result.json').write_text(json.dumps(result,indent=2)+'\n');manifest={'tag':tag,'xeff':args.xeff,'workflow_stage':'two_peak','selection':t['selection'],'threshold':thr,'threshold_source':t['thresholds_file'],'threshold_sha256':t['thresholds_sha256'],'inputs':records,'mass_range':[3.62,4.0],'mc_ranges':{k:v['range'] for k,v in mcresults.items()},'mean_half_range':.01,'scale_range':[.5,1.5],'scale_fixed':fixed_width,'scale_value':1.0 if fixed_width else None,'shape_strategy':result['shape_strategy'],'strategies':{'mc':2,'alt':2,'null':2},'null_definition':'Remove tested peak only; other peak and background remain profiled','input_contract_sha256':hashlib.sha256((art/'fit_input_contract.py').read_bytes()).hexdigest(),'significance_interpretation':'uncalibrated sqrt(q0), other peak and background profiled','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)

plot=(repo/'fitER/models/TwoPeakPlots.C').read_text()
if args.ppref:
 assert plot.count('PbPb DATA')==1;plot=plot.replace('PbPb DATA','ppRef DATA',1)
(art/'TwoPeakPlots.C').write_text(plot)
R.gSystem.SetBuildDir(str(art),True);assert R.gSystem.CompileMacro(str(art/'TwoPeakPlots.C'),'k')>0
R.PlotTwoData(str(out/'data_fit.pdf'),data,model,pdfs['x'],pdfs['psi'],bg,mass,alt,result['signals']['x']['Z_PL'])
for key,label in [('x','X'),('psi','#psi(2S)')]:
 info=mcresults[key];R.PlotPeakMC(str(out/(key+'_mc_template_fit.pdf')),mc[key],pdfs[key+'_mc'],mass,fitresults[key+'_mc'],info['mean'],info['sigma1'],info['sigma2'],info['fraction'],*info['range'],label)

if not args.ppref:subprocess.run([sys.executable,str(repo/'scripts/workflows/update_fit_reports.py')],check=True)

# A requested two-peak fit includes the conditional test after first-step quality checks.
subprocess.run([sys.executable,str(repo/'scripts/workflows/run_pb23_conditional_fit.py'),'--source',str(out)]+(['--no-report','--data-label','ppRef DATA'] if args.ppref else []),check=True)
