from pathlib import Path
import argparse,json,hashlib,shutil,subprocess,csv,re
import numpy as np,uproot
ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--resume',action='store_true');args=ap.parse_args();repo=Path.cwd();out=Path(args.output);out.mkdir(parents=True,exist_ok=args.resume);scratch=out/'artifacts';scratch.mkdir(exist_ok=args.resume)
if not args.resume:shutil.copy2(repo/'fitER/models/PbPbXEfficiencyFit.C',scratch)
assert (scratch/'PbPbXEfficiencyFit.C').read_bytes()==(repo/'fitER/models/PbPbXEfficiencyFit.C').read_bytes()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,a):p.write_text(json.dumps(a,indent=2)+'\n')
tags=['X_pb23_v7_fid5_10v2_rw0_xgb_v1']+['X_pb23_v3_fid3_'+v+'_rwr6range5v1_xgb_v1' for v in ['6v5','8v5','10v1']];labels=['10v2','6v5','8v5','10v1'];records=[];rows=[]
for tag,label in zip(tags,labels):
 base=repo.parent/'XGBoost/output/selected'/tag;tp=base/'cut_scan/weighted_signal_efficiency/thresholds.json';ts=json.loads(tp.read_text());points=[p for p in ts['thresholds'] if any(abs(p['target_efficiency']-e/100)<1e-8 for e in range(15,41,5))];assert len(points)==6;points.sort(key=lambda p:p['target_efficiency']);weight=ts['weight_branch'];weight=weight[0] if isinstance(weight,list) else weight
 rec={'tag':tag,'label':label,'thresholds_file':str(tp),'thresholds_sha256':sha(tp),'weight':weight,'points':points}
 if label!='10v2':
  mp=base/'fit_scan_manifest.pb23_pb24_simultaneous_mc_shape_nominal_v2.json';manifest=json.loads(mp.read_text());cat=manifest['pairing']['categories']['pb23'];cachebase=repo/'fitER/results/archive/earlier_workflows/pbpb_x_simultaneous_year_fit'/tag/'mc_shape_nominal_v2_fit_only_sqrtq0';ctx=json.loads((cachebase/'run_context.json').read_text());assert ctx['input_manifest_sha256']==sha(mp);assert cat['signal_mc']['event_weight_branch']==weight;assert cat['threshold_provenance']['sha256']==sha(tp)
  rec.update(manifest=str(mp),manifest_sha256=sha(mp),selection=cat['fiducial_selection']['expression']);files={k:cachebase/'cache/pb23'/v for k,v in [('data','DATA_fit_cache.root'),('mc','MC_fit_cache.root')]};trees={'data':cat['data']['tree'],'mc':cat['signal_mc']['tree']}
  for p in points:
   wp=next(q for q in manifest['working_points'] if abs(q['target_weighted_efficiency']-p['target_efficiency'])<1e-8);assert wp['categories']['pb23']['threshold']==p['score_threshold']
 elif args.resume:
  previous=json.loads((out/'scan_manifest.json').read_text());rec=next(t for t in previous['trainings'] if t['tag']==tag);assert rec['thresholds_sha256']==sha(tp);files={k:Path(v['path']) for k,v in rec['files'].items()};trees={k:v['tree'] for k,v in rec['files'].items()}
  for k,v in rec['files'].items():assert sha(files[k])==v['sha256']
 else:
  bs=base/'batch_apply_summary.json';b=json.loads(bs.read_text());assert b['input_datasets']['dataset_year']=='2023';rec.update(batch_apply_summary_sha256=sha(bs),selection=b['draw_selection']['fiducial_cut']['expression']);files={};trees={'data':'ntmix','mc':'ntmix_X3872'}
  # All scored entries were already fiducial-selected; verify selection exactly.
  for k,prefix in [('data','DATA'),('mc','MC')]:
   src=base/(prefix+'_with_score.root');stat=src.stat()
   with uproot.open(src) as f:
    a=f[trees[k]].arrays(['Bmass','Prediction','Bpt','By','BQvalue','Btrk2dR']+([weight] if k=='mc' else []),library='np')
   fid=(a['Bpt']>10)&(a['Bpt']<50)&(abs(a['By'])<1.6)&(a['BQvalue']<.15)&(a['Btrk2dR']>0)&(a['Btrk2dR']<.25);assert fid.all();mask=fid&(a['Bmass']>3.8)&(a['Bmass']<3.94)&(a['Prediction']>min(p['score_threshold'] for p in points));keep={n:a[n][mask] for n in ['Bmass','Prediction']+([weight] if k=='mc' else [])};keep['source_entry']=np.flatnonzero(mask).astype('int64');files[k]=scratch/(label+'_'+prefix+'_cache.root')
   with uproot.recreate(files[k]) as f:f.mktree(trees[k],{n:v.dtype for n,v in keep.items()});f[trees[k]].extend(keep)
   assert src.stat().st_mtime_ns==stat.st_mtime_ns;rec[k+'_source']={'path':str(src),'size':stat.st_size,'mtime_ns':stat.st_mtime_ns,'entries':len(mask),'cached':int(mask.sum())}
 rec['files']={k:{'path':str(p),'sha256':sha(p),'tree':trees[k]} for k,p in files.items()};records.append(rec);dump(out/'scan_manifest.json',{'schema_version':1,'contract':'pb23_four_training_independent_mc_shape_scan','macro_sha256':sha(scratch/'PbPbXEfficiencyFit.C'),'mass_range':[3.8,3.94],'mc_shape_range':[3.84,3.90],'width_scale_range':[.9,1.5],'trainings':records})
 for p in points:
  eff=round(100*p['target_efficiency']);point=out/label/f'xeff{eff}';point.mkdir(parents=True,exist_ok=args.resume);sel='Prediction > %.17g'%p['score_threshold'];strings=[label+f'_pb23_xeff{eff}',str(files['data']),trees['data'],str(files['mc']),trees['mc'],sel,sel,weight];expr='PbPbXEfficiencyFit.C+('+','.join(json.dumps(v) for v in strings)+',3.8,3.94,3.87169,0.005,0.9,1.5,28,'+json.dumps(str(point))+')'
  code=0
  if not (args.resume and (point/'fit_result.json').exists()):
   with (point/'fit.log').open('w') as log:code=subprocess.run(['root','-l','-b','-q',expr],cwd=scratch,stdout=log,stderr=subprocess.STDOUT).returncode
  assert code==0,(label,eff,code);result=json.loads((point/'fit_result.json').read_text());row={'tag':tag,'label':label,'xeff':eff,'threshold':p['score_threshold'],'weight_branch':weight,**result};rows.append(row);dump(out/'fit_summary.json',rows);print(label,eff,'Z',row['local_significance'],'yield',row['signal_yield'],'status',row['fit_status'],row['cov_qual'],flush=True)
 for k,info in rec['files'].items():assert sha(files[k])==info['sha256']
fields=list(rows[0]);
with (out/'fit_summary.csv').open('w') as f:w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows)
dump(out/'maxima.json',[max([r for r in rows if r['label']==label],key=lambda r:r['local_significance']) for label in labels]);print('ALL 24 COMPLETE',flush=True)
