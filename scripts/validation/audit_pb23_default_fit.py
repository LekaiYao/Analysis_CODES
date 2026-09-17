from pathlib import Path
import json,math,subprocess
import numpy as np,uproot,ROOT as R
import argparse
ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
R.gROOT.SetBatch(True);p=args.output;rows=json.loads((p/'fit_summary.json').read_text());manifest=json.loads((p/'scan_manifest.json').read_text());expected_points={(t['tag'],round(100*q['target_efficiency'])) for t in manifest['trainings'] for q in t['points']};assert len(rows)==len(expected_points);assert {(r['tag'],r['xeff']) for r in rows}==expected_points
checks=[]
for r in rows:
 info=next(t for t in manifest['trainings'] if t['tag']==r['tag']);d=p/r['label']/('xeff'+str(r['xeff']));f=R.TFile.Open(str(d/'fit_workspace.root'));q={k:f.Get('fit_result_'+k) for k in ['mc','alt','null']};assert all(q.values());z=math.sqrt(max(0,2*(q['null'].minNll()-q['alt'].minNll())));assert abs(z-r['local_significance'])<1e-8
 statuses={k:{'status':v.status(),'covQual':v.covQual(),'edm':v.edm()} for k,v in q.items()}
 expected={}
 for key in ['data','mc']:
  spec=info['files'][key]
  with uproot.open(spec['path']) as inp:a=inp[spec['tree']].arrays(['Bmass','Prediction']+([info['weight']] if key=='mc' else []),library='np')
  mask=(a['Prediction']>r['threshold'])&(a['Bmass']>manifest['mass_range'][0])&(a['Bmass']<manifest['mass_range'][1]);expected[key]=int(mask.sum());assert expected[key]==r['data_entries' if key=='data' else 'signal_mc_entries']
  if key=='mc':assert abs(a[info['weight']][mask].sum()-r['signal_mc_sumw'])<1e-6*max(1,r['signal_mc_sumw'])
 # exact Chebyshev positivity over the full mass support.
 c0,c1=r['chebyshev_a0'],r['chebyshev_a1'];grid=[-1.,1.]
 if c1 and -1 < -c0/(4*c1) < 1:grid.append(-c0/(4*c1))
 pos=min(1+c0*x+c1*(2*x*x-1) for x in grid);assert pos>0
 for pdf in d.glob('*.pdf'):assert 'Pages:' in subprocess.check_output(['pdfinfo',str(pdf)],text=True)
 valid=all(v['status']==0 and v['covQual']==3 and v['edm']<1e-3 for v in statuses.values());checks.append({'label':r['label'],'xeff':r['xeff'],'fits':statuses,'quality_pass':valid,'positive_background_min':pos,'q0_recomputed':z*z});f.Close()
result={'status':'PASS' if all(x['quality_pass'] for x in checks) else 'WARN_FIT_QUALITY','checks':checks,'count_weight_q0_pdf_checks':'PASS','toys':0,'bootstrap':0,'trials_correction':False};(p/'audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
