#!/usr/bin/env python3
"""Default PbPb23 single-year MC-shape efficiency scan; original fit macro unchanged."""
from pathlib import Path
import argparse,json,hashlib,shutil,subprocess,csv,re
import sys
from fit_input_contract import MASS_RANGE, POINTS, REVISION, load_training, selection_expression, verify_resume
ap=argparse.ArgumentParser()
ap.add_argument('--tag',action='append',required=True)
ap.add_argument('--output',help='Explicit single-tag staging directory; normally omit')
ap.add_argument('--variant',default='pb23_scan_wide',help='Semantic variant directory below the full ML tag')
ap.add_argument('--dry-run',action='store_true',help='Show destinations without reading data or fitting')
ap.add_argument('--resume',action='store_true')
ap.add_argument('--thresholds',type=Path,help='Validated downstream threshold supplement for explicitly requested extra points')
ap.add_argument('--xeff',nargs='+',type=int,default=POINTS,help='Explicit efficiency percentages; default remains 15 20 25 30 35 40')
ap.add_argument('--width-scale-range',nargs=2,type=float,default=[0.5,1.5],metavar=('MIN','MAX'))
ap.add_argument('--mean-half-range',type=float,default=0.01)
args=ap.parse_args();repo=Path(__file__).resolve().parents[2]
if len(set(args.xeff))!=len(args.xeff) or any(not 0<x<100 for x in args.xeff):ap.error('Efficiencies must be distinct integer percentages between 0 and 100')
scan_points=sorted(args.xeff)
if len(set(args.tag))!=len(args.tag):ap.error('Duplicate tags')
if not all(re.fullmatch(r'X_pb23_[A-Za-z0-9_]+',t) for t in args.tag):ap.error('Invalid PbPb23 X tag')
if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',args.variant):ap.error('Use a short semantic variant name without path separators')
if args.output and len(args.tag)!=1:ap.error('--output is only supported for one tag; omit it for automatic per-tag routing')
if not args.output:
 targets=[repo/'fitER/results/ml_fits'/tag/args.variant for tag in args.tag]
 for tag,target in zip(args.tag,targets):print(tag,'->',target,flush=True)
 if args.dry_run:print('stage=single_x; mass_range=',MASS_RANGE,'xeff=',scan_points,'metric=uncalibrated sqrt(q0)');sys.exit(0)
 for target in targets:
  if target.exists() and not args.resume:ap.error(str(target)+' already exists; choose a meaningful --variant or use --resume with unchanged settings')
 for tag,target in zip(args.tag,targets):
  cmd=[sys.executable,str(Path(__file__).resolve()),'--tag',tag,'--output',str(target),'--width-scale-range',*[str(x) for x in args.width_scale_range],'--mean-half-range',str(args.mean_half_range)]
  cmd.extend(['--xeff',*[str(x) for x in scan_points]])
  if args.thresholds:cmd.extend(['--thresholds',str(args.thresholds.resolve())])
  if args.resume:cmd.append('--resume')
  subprocess.run(cmd,cwd=repo,check=True)
 sys.exit(0)
out=Path(args.output).resolve()
results=repo/'fitER/results'
if out.is_relative_to(results) and (out.parent!=results/'ml_fits'/args.tag[0] or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*',out.name)):
 ap.error('Published results must use fitER/results/ml_fits/<full tag>/<variant>')
if args.dry_run:print(args.tag[0],'->',out);sys.exit(0)
import numpy as np,uproot,numexpr
assert 0 < args.width_scale_range[0] < args.width_scale_range[1]
assert 0 < args.mean_half_range < 0.03
prior=json.loads((out/'scan_manifest.json').read_text()) if args.resume else None
if prior is not None:
 assert prior.get('workflow_revision')==REVISION and prior['mass_range']==MASS_RANGE, 'Old mass range/workflow cannot resume; choose a new variant'
out.mkdir(parents=True,exist_ok=args.resume);scratch=out/'artifacts';scratch.mkdir(exist_ok=args.resume)
if not args.resume:
 shutil.copy2(repo/'fitER/models/PbPbXEfficiencyFit.C',scratch)
 shutil.copy2(__file__,scratch)
 shutil.copy2(repo/'scripts/workflows/fit_input_contract.py',scratch)
else:
 for name in ['run_pb23_default_fit.py','fit_input_contract.py']:
  assert (scratch/name).read_bytes()==(repo/'scripts/workflows'/name).read_bytes(), 'Resume code changed'

assert (scratch/'PbPbXEfficiencyFit.C').read_bytes()==(repo/'fitER/models/PbPbXEfficiencyFit.C').read_bytes()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,a):p.write_text(json.dumps(a,indent=2)+'\n')
tags=args.tag;labels=[t.split('_fid',1)[1].split('_')[1] for t in tags];assert len(set(labels))==len(labels);assert all(t.startswith('X_pb23_') for t in tags);records=[];rows=[]
for tag,label in zip(tags,labels):
 rec,allpoints=load_training(repo,tag);rec['label']=label
 if args.thresholds:
  extra=json.loads(args.thresholds.read_text());assert extra['schema_version']==1 and extra['train_tag']==tag
  assert extra['parent_thresholds_sha256']==rec['thresholds_sha256'] and extra['selection']==rec['selection']
  assert extra['weight_branch']==[rec['weight']] and extra['comparison']=='Prediction > threshold'
  assert extra['existing_thresholds_exactly_reproduced']==len(allpoints)
  assert sha(Path(extra['reference_file']))==extra['reference_sha256']
  rec['parent_thresholds_file']=rec['thresholds_file'];rec['parent_thresholds_sha256']=rec['thresholds_sha256']
  rec['thresholds_file']=str(args.thresholds.resolve());rec['thresholds_sha256']=sha(args.thresholds)
  allpoints=extra['thresholds']

 points=[p for p in allpoints if any(abs(p['target_efficiency']-e/100)<1e-8 for e in scan_points)];assert len(points)==len(scan_points);points.sort(key=lambda p:p['target_efficiency']);weight=rec['weight'];rec['points']=points
 if args.resume:
  rec=verify_resume(prior,rec,points,args.width_scale_range,args.mean_half_range);files={k:Path(v['path']) for k,v in rec['files'].items()};trees={k:v['tree'] for k,v in rec['files'].items()}
 else:
  files={};trees={k:v['tree'] for k,v in rec['sources'].items()}
  for k,prefix in [('data','DATA'),('mc','MC')]:
   src=Path(rec['sources'][k]['path']);stat=src.stat()
   with uproot.open(src) as f:
    a=f[trees[k]].arrays(sorted(set(['Bmass','Prediction']+re.findall(r'\bB\w+\b',rec['selection'])+([weight] if k=='mc' else []))),library='np')
   fid=numexpr.evaluate(selection_expression(rec['selection']),local_dict=a);assert fid.all(), 'Scored input does not satisfy declared fiducial selection'
   mask=fid&(a['Bmass']>MASS_RANGE[0])&(a['Bmass']<MASS_RANGE[1])&(a['Prediction']>min(p['score_threshold'] for p in points));keep={n:a[n][mask] for n in ['Bmass','Prediction']+([weight] if k=='mc' else [])};keep['source_entry']=np.flatnonzero(mask).astype('int64');files[k]=scratch/(label+'_'+prefix+'_cache.root')
   with uproot.recreate(files[k]) as f:f.mktree(trees[k],{n:v.dtype for n,v in keep.items()});f[trees[k]].extend(keep)
   assert src.stat().st_mtime_ns==stat.st_mtime_ns;rec[k+'_source']={'path':str(src),'size':stat.st_size,'mtime_ns':stat.st_mtime_ns,'entries':len(mask),'cached':int(mask.sum())}
 rec['files']={k:{'path':str(p),'sha256':sha(p),'tree':trees[k]} for k,p in files.items()};records.append(rec);dump(out/'scan_manifest.json',{'schema_version':2,'workflow_revision':REVISION,'strategies':{'mc':2,'alt':2,'null':2},'significance_interpretation':'uncalibrated sqrt(q0), background profiled','contract':'pb23_default_independent_mc_shape_scan','macro_sha256':sha(scratch/'PbPbXEfficiencyFit.C'),'target_efficiencies_percent':scan_points,'mass_range':MASS_RANGE,'mc_shape_range':[3.84,3.90],'width_scale_range':args.width_scale_range,'mean_nominal':3.87169,'mean_half_range':args.mean_half_range,'mean_range':[3.87169-args.mean_half_range,3.87169+args.mean_half_range],'trainings':records})
 for p in points:
  eff=round(100*p['target_efficiency']);point=out/label/f'xeff{eff}';point.mkdir(parents=True,exist_ok=args.resume);sel='Prediction > %.17g'%p['score_threshold'];strings=[label+f'_pb23_xeff{eff}',str(files['data']),trees['data'],str(files['mc']),trees['mc'],sel,sel,weight];expr='PbPbXEfficiencyFit.C+('+','.join(json.dumps(v) for v in strings)+(',3.75,4.00,3.87169,%.17g,%.17g,%.17g,50,'%(args.mean_half_range,*args.width_scale_range))+json.dumps(str(point))+')'
  code=0
  if not (args.resume and (point/'fit_result.json').exists()):
   with (point/'fit.log').open('w') as log:code=subprocess.run(['root','-l','-b','-q',expr],cwd=scratch,stdout=log,stderr=subprocess.STDOUT).returncode
  assert code==0,(label,eff,code);result=json.loads((point/'fit_result.json').read_text());row={'tag':tag,'label':label,'xeff':eff,'threshold':p['score_threshold'],'weight_branch':weight,**result};rows.append(row);dump(out/'fit_summary.json',rows);print(label,eff,'Z',row['local_significance'],'yield',row['signal_yield'],'status',row['fit_status'],row['cov_qual'],flush=True)
 for k,info in rec['files'].items():assert sha(files[k])==info['sha256']
fields=list(rows[0]);
with (out/'fit_summary.csv').open('w') as f:w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows)
dump(out/'maxima.json',[max([r for r in rows if r['label']==label],key=lambda r:r['local_significance']) for label in labels]);print('ALL',len(rows),'COMPLETE',flush=True)

subprocess.run([sys.executable,str(repo/'scripts/validation/audit_pb23_default_fit.py'),'--output',str(out)],check=True)

pdfs=[str(out/row['label']/f"xeff{row['xeff']}"/'data_fit.pdf') for row in rows]
subprocess.run(['pdfunite',*pdfs,str(out/'all_fits.pdf')],check=True)
subprocess.run([sys.executable,str(repo/'scripts/workflows/update_fit_reports.py')],check=True)
