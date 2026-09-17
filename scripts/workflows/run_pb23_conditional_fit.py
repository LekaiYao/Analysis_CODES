from pathlib import Path
import json,math,hashlib,shutil,argparse,subprocess,sys
ap=argparse.ArgumentParser(description='Explicit conditional X test from a completed two-peak fit')
ap.add_argument('--source',type=Path,required=True,help='Directory of the selected two-peak result');ap.add_argument('--dry-run',action='store_true');ap.add_argument('--no-report',action='store_true');ap.add_argument('--data-label',choices=['PbPb DATA','ppRef DATA'],default='PbPb DATA');args=ap.parse_args()
repo=Path(__file__).resolve().parents[2];base=args.source.resolve();out=base/'conditional'
first_result=json.loads((base/'fit_result.json').read_text());assert first_result['mass_range']==[3.62,4.0]
assert 'x' in first_result['signals'] and 'psi' in first_result['signals']
if args.dry_run:print(out,'reuse step1; alt floats X/background yields; null floats background yield');raise SystemExit(0)
assert first_result['audit']['quality_pass'], 'Review failed first-step fit before using conditional test'
assert not out.exists(), 'Refusing to overwrite existing result'
import ROOT as R
R.gROOT.SetBatch(True);R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
out.mkdir();art=out/'artifacts';art.mkdir();shutil.copy2(__file__,art)
src=base/'fit_workspace.root';f=R.TFile.Open(str(src));w=f.Get('ws_nominal');first=f.Get('fit_result_alt');w.allVars().assignValueOnly(first.floatParsFinal());model=w.pdf('model');data=w.data('data_data');mass=w.var('Bmass');params=model.getParameters(data)
for p in params:p.setConstant(True)
for n in ['x_yield','nbkg']:w.var(n).setConstant(False)
fixed={p.GetName():p.getVal() for p in params if p.isConstant()};initial={p.GetName():p.getVal() for p in params}
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
def fit(pdf):return fit_with_one_restart(pdf,data,[R.RooFit.Save(),R.RooFit.Extended(True),R.RooFit.Range('all'),R.RooFit.Strategy(2),R.RooFit.Hesse(True),R.RooFit.PrintLevel(-1),R.RooFit.Warnings(False),R.RooFit.Verbose(False)])
def info(r):return {'nll':r.minNll(),'status':r.status(),'covQual':r.covQual(),'edm':r.edm(),'floating':{p.GetName():{'value':p.getVal(),'error':p.getError(),'range':[p.getMin(),p.getMax()]} for p in r.floatParsFinal()}}
alt=fit(model);assert {p.GetName() for p in alt.floatParsFinal()}=={'x_yield','nbkg'};altinfo=info(alt);snap=params.snapshot()
null=R.RooAddPdf('model_null_x','psi2S plus background; X removed',R.RooArgList(w.pdf('psi_signal'),w.pdf('backgroundPdf')),R.RooArgList(w.var('psi_yield'),w.var('nbkg')));nullfit=fit(null);assert {p.GetName() for p in nullfit.floatParsFinal()}=={'nbkg'};nullinfo=info(nullfit)
for n,val in fixed.items():assert w.var(n).getVal()==val
nll0=null.createNLL(data,R.RooFit.Extended(True),R.RooFit.Range('all'));assert abs(nll0.getVal()-nullfit.minNll())<1e-6
w.var('x_yield').setVal(0);nll=model.createNLL(data,R.RooFit.Extended(True),R.RooFit.Range('all'));assert abs(nll.getVal()-nll0.getVal())<1e-6
params.assignValueOnly(snap);assert abs(nll.getVal()-alt.minNll())<1e-6;q0=max(0.,2*(nullfit.minNll()-alt.minNll()));z=math.sqrt(q0)
r={'tag':first_result['tag'],'xeff':first_result['xeff'],'interpretation':'conditional same-data fixed-parameter sqrt(q0); calibration not established','mass_range':[3.62,4.0],'data_entries':data.numEntries(),'step1_reused':{'path':str(src),'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'fit':'fit_result_alt','nll':first.minNll()},'fixed_parameters':fixed,'initial_parameters':initial,'alt':altinfo,'null':nullinfo,'q0':q0,'Z_PL':z,'audit':{'quality_pass':all(x['status']==0 and x['covQual']==3 and x['edm']<1e-3 for x in [altinfo,nullinfo]),'fixed_parameters_unchanged':True,'null_removed_x_matches_zero_yield':True,'nll_recomputed':True}}
output=R.TFile.Open(str(out/'fit_workspace.root'),'RECREATE');ww=R.RooWorkspace('ws_conditional');imp=getattr(ww,'import');imp(data);imp(model);imp(null,R.RooFit.RecycleConflictNodes());ww.Write();alt.Write('fit_result_alt');nullfit.Write('fit_result_null');first.Write('fit_result_step1');output.Close();(out/'fit_result.json').write_text(json.dumps(r,indent=2)+'\n')
plot=(repo/'fitER/models/TwoPeakPlots.C').read_text();assert plot.count('PbPb DATA')==1;plot=plot.replace('PbPb DATA',args.data_label,1);a='fit.floatParsFinal().find("psi_yield")';assert plot.count(a)==1;plot=plot.replace(a,'plotParameters->find("psi_yield")',1);a='"N_{#psi}=%.1f #pm %.1f",psiY->getVal(),psiY->getError()';assert plot.count(a)==1;plot=plot.replace(a,'"N_{#psi}=%.1f (fixed)",psiY->getVal()',1);(art/'TwoPeakPlots.C').write_text(plot)
R.gSystem.SetBuildDir(str(art),True);assert R.gSystem.CompileMacro(str(art/'TwoPeakPlots.C'),'k')>0;R.PlotTwoData(str(out/'data_fit.pdf'),data,model,w.pdf('x_signal'),w.pdf('psi_signal'),w.pdf('backgroundPdf'),mass,alt,z)
print(json.dumps(r,indent=2),flush=True)

if not args.no_report:subprocess.run([sys.executable,str(repo/'scripts/workflows/update_fit_reports.py')],check=True)
