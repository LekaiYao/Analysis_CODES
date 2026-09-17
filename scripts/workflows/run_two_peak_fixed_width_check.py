"""First-step-only two-peak check: float means, fix both MC width scales to 1."""
from pathlib import Path
import argparse,json,hashlib,math,shutil
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',type=Path,required=True);ap.add_argument('--dry-run',action='store_true');args=ap.parse_args()
repo=Path(__file__).resolve().parents[2];base=args.source.resolve();old=json.loads((base/'fit_result.json').read_text());manifest=json.loads((base/'manifest.json').read_text())
assert old['mass_range']==[3.62,4.0] and old['audit']['quality_pass']
assert base.parent.name=='two_peak' and base.name==f"xeff{old['xeff']}"
out=base.parents[1]/'checks/two_peak_mc_width'/base.name
if args.dry_run:print(out);raise SystemExit(0)
assert not out.exists(),'Refusing to overwrite existing result'
import ROOT as R
R.gROOT.SetBatch(True);R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
source=base/'fit_workspace.root';source_hash=hashlib.sha256(source.read_bytes()).hexdigest();f=R.TFile.Open(str(source));w=f.Get('ws_nominal');data=w.data('data_data');model=w.pdf('model');mass=w.var('Bmass');n=data.numEntries();assert n==old['data_entries']
# Reset the same original initial values; the only physical change is fixed scale=1.
for key,nom in [('x',3.87169),('psi',3.68610)]:
 for name in ['sigma1','sigma2','fraction']:
  v=w.var(key+'_'+name);assert v.isConstant();assert abs(v.getVal()-old['mc'][key][name])<1e-12
 mean=w.var(key+'_mean');assert abs(mean.getMin()-(nom-.01))<1e-10 and abs(mean.getMax()-(nom+.01))<1e-10
 mean.setVal(old['mc'][key]['mean']);mean.setConstant(False)
 scale=w.var(key+'_scale');scale.setVal(1.);scale.setConstant(True)
 near=data.sumEntries(f'abs(Bmass-{nom})<0.005');y=w.var(key+'_yield');y.setVal(max(0.,.4*near));y.setConstant(False)
for name,value in [('a0',-.35),('a1',-.05),('nbkg',.7*n)]:w.var(name).setVal(value);w.var(name).setConstant(False)
params=model.getParameters(data);initial={p.GetName():{'value':p.getVal(),'fixed':p.isConstant()} for p in params}
fixed={p.GetName():p.getVal() for p in params if p.isConstant()}
out.mkdir(parents=True);art=out/'artifacts';art.mkdir();shutil.copy2(__file__,art)
def record(q):return {'status':q.status(),'covQual':q.covQual(),'edm':q.edm(),'nll':q.minNll(),'history':[(q.statusLabelHistory(i),q.statusCodeHistory(i)) for i in range(q.numStatusHistory())],'parameters':{p.GetName():{'value':p.getVal(),'error':p.getError()} for p in q.floatParsFinal()}}
def passed(q):return q.status()==0 and q.covQual()==3 and math.isfinite(q.edm()) and q.edm()<.001
opts=[R.RooFit.Save(),R.RooFit.Extended(True),R.RooFit.Range('all'),R.RooFit.Strategy(2),R.RooFit.Hesse(True),R.RooFit.PrintLevel(-1),R.RooFit.Warnings(False),R.RooFit.Verbose(False)]
fit=model.fitTo(data,*opts);attempts=[record(fit)];first=None
if math.isfinite(fit.minNll()) and not passed(fit):
 first=fit;params.assignValueOnly(fit.floatParsFinal());fit=model.fitTo(data,*opts);attempts.append(record(fit))
assert {p.GetName() for p in fit.floatParsFinal()}=={'x_mean','x_yield','psi_mean','psi_yield','nbkg','a0','a1'}
for name,value in fixed.items():assert w.var(name).isConstant() and w.var(name).getVal()==value
nll=model.createNLL(data,R.RooFit.Extended(True),R.RooFit.Range('all'));assert abs(nll.getVal()-fit.minNll())<1e-6
aa=w.var('a0').getVal();bb=w.var('a1').getVal();xx=[-1.,1.]
if bb and -1<-aa/(4*bb)<1:xx.append(-aa/(4*bb))
minimum=min(1+aa*x+bb*(2*x*x-1) for x in xx);assert minimum>0
assert hashlib.sha256(source.read_bytes()).hexdigest()==source_hash
for inp in manifest['inputs'].values():assert hashlib.sha256(Path(inp['cache']).read_bytes()).hexdigest()==inp['sha256']
def param(v):return {'value':v.getVal(),'error':None if v.isConstant() else v.getError(),'fixed':v.isConstant(),'range':[v.getMin(),v.getMax()],'at_boundary':not v.isConstant() and min(v.getVal()-v.getMin(),v.getMax()-v.getVal())<=1e-4*(v.getMax()-v.getMin())}
result={'tag':old['tag'],'xeff':old['xeff'],'mass_range':old['mass_range'],'data_entries':n,'mc':old['mc'],'signals':{k:{name:param(w.var(k+'_'+name)) for name in ['mean','scale','yield']} for k in ['x','psi']},'background':{name:param(w.var(name)) for name in ['a0','a1','nbkg']},'fits':{'alt':record(fit)},'audit':{'quality_pass':passed(fit),'fixed_parameters_unchanged':True,'nll_recomputed':True,'source_hash_unchanged':True,'cache_hashes_pass':True,'background_min':minimum},'first_step_only':True,'conditional_test_performed':False,'comparison':{'original_float_shape_nll':old['fits']['alt']['nll'],'fixed_width_minus_float_nll':fit.minNll()-old['fits']['alt']['nll']}}
(out/'fit_result.json').write_text(json.dumps(result,indent=2)+'\n');(out/'retry_history.json').write_text(json.dumps(attempts,indent=2)+'\n')
meta={'source_workspace':str(source),'source_sha256':source_hash,'source_manifest_sha256':hashlib.sha256((base/'manifest.json').read_bytes()).hexdigest(),'selection':manifest['selection'],'threshold':manifest['threshold'],'mc_reused':True,'mc_width_source':'original same-tag xeff50 templates, not common-precut templates','fixed_parameters':fixed,'initial_parameters':initial,'scope':'alt first step only; no null or conditional fits','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(out/'manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
output=R.TFile.Open(str(out/'fit_workspace.root'),'RECREATE');w.Write();fit.Write('fit_result_alt')
for k in ['x','psi']:f.Get('fit_result_'+k+'_mc').Write('fit_result_'+k+'_mc')
if first:first.Write('fit_result_alt_first_attempt')
output.Close()
plot=(repo/'fitER/models/TwoPeakPlots.C').read_text();line='    stats.AddText(Form("Z_{PL,X}=%.3f", z));';assert plot.count(line)==1;plot=plot.replace(line,'    stats.AddText("MC widths fixed; mean free");',1)
(art/'TwoPeakPlots.C').write_text(plot);R.gSystem.SetBuildDir(str(art),True);assert R.gSystem.CompileMacro(str(art/'TwoPeakPlots.C'),'k')>0
R.PlotTwoData(str(out/'data_fit.pdf'),data,model,w.pdf('x_signal'),w.pdf('psi_signal'),w.pdf('backgroundPdf'),mass,fit,0.)
print(json.dumps(result,indent=2),flush=True)
