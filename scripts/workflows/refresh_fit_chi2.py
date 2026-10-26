"""Explicitly authorized in-place combined reruns / individual DATA plot refresh.
Run in ROOT 6.32.02. Back up every replaced artifact; do not update slides requests.
"""
from pathlib import Path
import argparse, json, hashlib, shutil, subprocess, sys, tempfile, os, math, datetime
import numpy as np
import uproot
ROOT_DIR=Path(__file__).resolve().parents[2]
AUDIT_DIR=ROOT_DIR/'fitER/results/comparisons/fit_chi2'
def read(p): return json.loads(Path(p).read_text())
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d): Path(p).write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def quality_audit(p):
 d=read(p);bins=d['bins'];values=[(b['observed']-b['expected_curve_average'])/(b['error_low'] if b['observed']>b['expected_curve_average'] else b['error_high']) for b in bins if b['observed']!=0]
 assert len(values)==d['participating_bins']
 assert math.isclose(sum(x*x for x in values),d['chi2'],rel_tol=1e-12,abs_tol=1e-10)
 assert d['ndf']==d['participating_bins']-d['floating_parameters']
 assert math.isclose(d['chi2']/d['ndf'],d['chi2_ndf'],rel_tol=1e-12)
 return {k:v for k,v in d.items() if k!='bins'}
def publish(stage,out,originals):
 records=[]
 for p in sorted(stage.rglob('*')):
  if not p.is_file():continue
  if 'build' in p.relative_to(stage).parts:continue
  target=out/p.relative_to(stage)
  if target.exists():
   assert sha(target)==originals[str(target)],str(target)+' changed during run'
   backup=AUDIT_DIR/'metadata/originals'/target.relative_to(ROOT_DIR)
   assert not backup.exists(),str(backup)+' already backed up; refusing duplicate refresh'
   backup.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(target,backup)
  else:backup=None
  target.parent.mkdir(parents=True,exist_ok=True)
  tmp=target.with_name(target.name+'.chi2_tmp');shutil.copy2(p,tmp);os.replace(tmp,target)
  assert sha(target)==sha(p)
  records.append({'path':str(target),'old_sha256':originals.get(str(target)),'new_sha256':sha(target),'backup':str(backup) if backup else None})
 return records

def run(mode,out):
 import ROOT as R
 R.gROOT.SetBatch(True);R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
 assert R.gROOT.GetVersion()=='6.32.02'
 out=out.resolve();root=ROOT_DIR/'fitER/results'/('combined_fits' if mode=='combined' else 'ml_fits');assert out.is_relative_to(root)
 originals={str(p):sha(p) for p in out.rglob('*') if p.is_file() and not any(x in p.parts for x in ['build','initial_local_minimum'])}
 stage=Path(tempfile.mkdtemp(prefix='fit_chi2_',dir='/tmp/leyao'));art=stage/'artifacts';art.mkdir()
 started=datetime.datetime.now(datetime.timezone.utc).isoformat()
 if mode=='combined':
  from combined_two_peak_model import CombinedModel
  from run_pb18_pb23_combined_fit import plots
  old=read(out/'fit_result.json');manifest=read(out/'manifest.json');arrays={}
  for y in ('2018','2023'):
   arrays[y]={}
   for k,spec in old['inputs'][y].items():
    p=out/spec['cache'];assert sha(p)==spec['cache_sha256']
    st=Path(spec['path']).stat();assert (st.st_size,st.st_mtime_ns)==(spec['size'],spec['mtime_ns'])
    with uproot.open(p) as f:arrays[y][k]=f[spec['tree']].arrays(library='np')
   for p,h in manifest['categories'][y]['provenance'].items():assert sha(p)==h,p
  model=CombinedModel(arrays);model.fit_mc();model.build_data_model()
  # Reuse the already accepted full optimum only as an initial point. All fits run again.
  for name,v in old['fits']['full']['floating'].items():model.params.find(name).setVal(v['value'])
  model.fit_full();model.fit_conditional()
  assert abs(model.full_nll['sum_year_nll']-old['full']['nll_audit']['sum_year_nll'])<.01
  model.save_workspace(stage/'fit_workspace.root');diagnostics=plots(model,stage,ROOT_DIR)
  for p in stage.glob('*.chi2.json'):quality_audit(p)
  report=model.report();report.update(status='PASS',workflow_revision=old['workflow_revision'],inputs=old['inputs'],data_entries=old['data_entries'],data_chi2=diagnostics)
  refresh={'started_UTC':started,'operation':'MC/full/H1/H0 refitted; previous accepted full parameters as initial values','input_cache_reused':True,'old_result_sha256':originals[str(out/'fit_result.json')],'old_significance':old['conditional']['significance'],'new_significance':report['conditional']['significance'],'delta_full_nll':model.full_nll['sum_year_nll']-old['full']['nll_audit']['sum_year_nll'],'chi2_definition':'ROOT 6.32.02 RooCurve nonzero bins, asymmetric Poisson errors, curve bin average; extended expected yield; combined subtracts 14 full / 2 conditional parameters once'}
  report['chi2_refresh']=refresh
  for n in ('combined_two_peak_model.py','run_pb18_pb23_combined_fit.py','refresh_fit_chi2.py'):shutil.copy2(ROOT_DIR/'scripts/workflows'/n,art/n)
  manifest['implementation'].update({p.name:sha(p) for p in art.iterdir() if p.suffix in ('.py','.C')})
  manifest['chi2_refresh']=refresh;manifest.pop('numerical_restart',None)
  save(stage/'fit_result.json',report);save(stage/'manifest.json',manifest)
  txt=(out/'README.md').read_text();txt=txt.replace(f"Combined significance: {old['conditional']['significance']:.6f}.",f"Combined significance: {report['conditional']['significance']:.6f}.")
  txt+='\n2026-09-30: MC/full/H1/H0 rerun with unchanged input/model; DATA χ², ndf and χ²/ndf in plots and data_chi2.json. Previous result backed up under comparisons/fit_chi2/metadata/originals.\n'
  (stage/'README.md').write_text(txt)
  numerical={'before':old['conditional']['significance'],'after':report['conditional']['significance'],'data_chi2':diagnostics}
 else:
  first=read(out/'fit_result.json');assert first['audit']['quality_pass']
  numerical={};files=[]
  for kind,sub,wsname in [('full','', 'ws_nominal'),('conditional','conditional','ws_conditional')]:
   src=out/sub;dest=stage/sub;dest.mkdir(exist_ok=True);(dest/'artifacts').mkdir(exist_ok=True)
   result=read(src/'fit_result.json');file=R.TFile.Open(str(src/'fit_workspace.root'));w=file.Get(wsname);assert w
   fit=file.Get('fit_result_alt');assert fit and fit.status()==0 and fit.covQual()==3
   w.allVars().assignValueOnly(fit.floatParsFinal());data=w.data('data_data');model=w.pdf('model');mass=w.var('Bmass')
   z=first['signals']['x']['Z_PL'] if kind=='full' else result['Z_PL']
   plot=(ROOT_DIR/'fitER/models/TwoPeakPlots.C').read_text();assert plot.count('PbPb DATA')==1
   year=str(first.get('year') or ('2018' if 'pb18' in first['tag'] else '2023'));plot=plot.replace('PbPb DATA','PbPb '+year+' DATA',1)
   macro=dest/'artifacts/TwoPeakPlots.C';macro.write_text(plot)
   build=Path('/tmp/leyao/fit_chi2/plot_cache')/hashlib.sha256(plot.encode()).hexdigest();build.mkdir(parents=True,exist_ok=True)
   cp=build/'TwoPeakPlots.C'
   if not cp.exists():cp.write_text(plot)
   assert cp.read_text()==plot;R.gSystem.SetBuildDir(str(build),True)
   if kind=='full':
    import fcntl
    with (build/'compile.lock').open('a') as lock:
     fcntl.flock(lock,fcntl.LOCK_EX)
     assert R.gSystem.CompileMacro(str(cp),'k')>0
   R.PlotTwoData(str(dest/'data_fit.pdf'),data,model,w.pdf('x_signal'),w.pdf('psi_signal'),w.pdf('backgroundPdf'),mass,fit,z)
   numerical[kind]=quality_audit(dest/'data_fit.chi2.json')
   save(dest/'artifacts/plot_build.json',{'cache':str(build),'root_version':R.gROOT.GetVersion(),'source_sha256':sha(macro),'plot_only':True})
   save(dest/'plot_refresh.json',{'started_UTC':started,'plot_only':True,'workspace_sha256':sha(src/'fit_workspace.root'),'fit_result_sha256':sha(src/'fit_result.json'),'source_sha256':sha(macro),'data_chi2':numerical[kind]})
   files.append(file)
  save(stage/'data_chi2.json',numerical)
  # No individual ROOT/fit JSON is written or refitted.
 changed=publish(stage,out,originals)
 audit={'status':'PASS','mode':mode,'directory':str(out),'started_UTC':started,'ended_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'numerical':numerical,'changed':changed}
 if mode=='individual':
  for sub in ('','conditional'):
   for n in ('fit_workspace.root','fit_result.json'):assert sha(out/sub/n)==originals[str(out/sub/n)]
  audit['individual_fit_and_workspace_unchanged']=True
 name=hashlib.sha256(str(out).encode()).hexdigest()[:16]+'.json';dest=AUDIT_DIR/'metadata/runs';dest.mkdir(parents=True,exist_ok=True);save(dest/name,audit)
 shutil.rmtree(stage)
 print(json.dumps({'status':'PASS','mode':mode,'directory':str(out),'audit':str(dest/name)}),flush=True)
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--mode',choices=['combined','individual'],required=True);ap.add_argument('--directory',type=Path,required=True);args=ap.parse_args();run(args.mode,args.directory)
