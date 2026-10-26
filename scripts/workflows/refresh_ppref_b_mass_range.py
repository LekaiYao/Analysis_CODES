#!/usr/bin/env python3
"""Stage all current ppRef B mass spectra at strict 5 < Bmass < 6 GeV."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,shutil,subprocess,re,copy
import numpy as np
import uproot
import ROOT as R
from ppref_b_plot_style import plot
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0);R.gStyle.SetTextFont(42)
assert R.gROOT.GetVersion()=='6.32.02'
repo=Path(__file__).resolve().parents[2];out=repo/'plotER/results/Bmass_ppRef';stage=Path('/tmp/leyao/ppref_b_mass_5to6');assert not stage.exists();stage.mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
oldhashes={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}
(stage/'before_hashes.json').write_text(json.dumps(oldhashes,indent=2)+'\n')
shutil.copytree(out,stage/'before')
pro=json.loads((out/'metadata/provenance.json').read_text());cutpro=json.loads((out/'metadata/2D4/provenance.json').read_text());cfgs=json.loads((out/'metadata/variables.json').read_text())
cuts=json.loads((out/'metadata/cuts.json').read_text());assert set(cuts)=={'baseline','2D4'}
for ch in pro['channels']:assert sorted(p.name for p in (out/ch/'pdf').glob('*.pdf'))==['2D4.pdf','baseline.pdf']
files={key:R.TFile(str(stage/(key+'_histograms.root')),'RECREATE') for key in cuts}
report={'status':'pending','created_UTC':datetime.now(timezone.utc).isoformat(),'display_range':[5.,6.],'display_expression':'Bmass > 5 && Bmass < 6','bins':100,'bin_width_GeV':.01,'selection_unchanged':cuts,'channels':{},'source_SHA_recomputed':False,'script_sha256':sha(Path(__file__)),'style_sha256':sha(Path(__file__).with_name('ppref_b_plot_style.py'))}
for ch,c in pro['channels'].items():
 src=Path(c['source']);before=src.stat();assert (before.st_size,before.st_mtime_ns)==(c['fingerprint']['size'],c['fingerprint']['mtime_ns'])
 stats={k:{'total':0,'finite':0,'nonfinite':0,'min':None,'max':None,'underflow':0,'overflow':0,'equal_lower_edge':0,'equal_upper_edge':0} for k in cuts};counts={k:np.zeros(100,dtype=np.int64) for k in cuts}
 with uproot.open(src) as f:
  for a in f[c['tree']].iterate(['Bmass','Bnorm_svpvDistance_2D'],step_size='80 MB',library='np'):
   assert np.isfinite(a['Bnorm_svpvDistance_2D']).all()
   for key,mask in [('baseline',np.ones(len(a['Bmass']),dtype=bool)),('2D4',a['Bnorm_svpvDistance_2D']>4)]:
    x=a['Bmass'][mask].astype(float);s=stats[key];s['total']+=len(x);ok=np.isfinite(x);s['nonfinite']+=int((~ok).sum());x=x[ok];s['finite']+=len(x)
    if len(x):s['min']=float(x.min()) if s['min'] is None else min(s['min'],float(x.min()));s['max']=float(x.max()) if s['max'] is None else max(s['max'],float(x.max()))
    inside=(x>5)&(x<6);counts[key]+=np.histogram(x[inside],np.linspace(5,6,101))[0];s['underflow']+=int((x<=5).sum());s['overflow']+=int((x>=6).sum());s['equal_lower_edge']+=int((x==5).sum());s['equal_upper_edge']+=int((x==6).sum())
 df=R.RDataFrame(c['tree'],str(src));actions={}
 for key,expr in [('baseline','true'),('2D4','Bnorm_svpvDistance_2D > 4')]:
  f=df.Filter(expr);inside=f.Filter('Bmass > 5 && Bmass < 6');actions[key]=(f.Count(),inside.Histo1D((ch+'_'+key+'_audit','',100,5.,6.),'Bmass'),f.Filter('Bmass <= 5').Count(),f.Filter('Bmass >= 6').Count())
 report['channels'][ch]={}
 for key in cuts:
  n,hp,low,high=actions[key];h=hp.GetValue();s=stats[key];assert int(n.GetValue())==s['total'];assert np.array_equal(counts[key],[h.GetBinContent(i) for i in range(1,101)]);assert np.allclose([h.GetBinError(i)**2 for i in range(1,101)],counts[key]);assert low.GetValue()==s['underflow'] and high.GetValue()==s['overflow'];assert sum(counts[key])+s['underflow']+s['overflow']==s['finite']
  expected=c['entries'] if key=='baseline' else cutpro['channels'][ch]['statistics']['total'];assert s['total']==expected
  cfg=copy.deepcopy(cfgs[ch]['Bmass']);cfg.update(range=[5.,6.],bins=100,range_source='user requested common range 5 < Bmass < 6',display_expression='Bmass > 5 && Bmass < 6')
  label='Raw input; no additional cuts' if key=='baseline' else 'Bnorm_svpvDistance_2D > 4'
  name='Bmass' if key=='baseline' else 'Bmass_2D4'
  plot(ch,name,cfg,counts[key],s,stage/f'{ch}_{key}.pdf',{'Bmass_ppRef':files[key]},mass=True,selection_label=label)
  report['channels'][ch][key]={'statistics':s,'bin_counts':counts[key].tolist(),'config':cfg,'in_range':int(sum(counts[key])),'ROOT_numpy_selection_bins_errors':'PASS'}
 subprocess.run(['pdfunite',str(stage/f'{ch}_baseline.pdf'),str(stage/f'{ch}_2D4.pdf'),str(stage/f'{ch}_all_cuts.pdf')],check=True)
 after=src.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
 print(ch,{k:int(sum(counts[k])) for k in cuts},flush=True)
for f in files.values():f.Close()
for key in cuts:
 with uproot.open(stage/(key+'_histograms.root')) as f:
  for ch in report['channels']:
   name=ch+('_Bmass' if key=='baseline' else '_Bmass_2D4');assert np.array_equal(f[name].values(),report['channels'][ch][key]['bin_counts']);assert np.allclose(f[name].variances(),report['channels'][ch][key]['bin_counts'])
for ch in report['channels']:
 for key in ['baseline','2D4','all_cuts']:
  p=stage/f'{ch}_{key}.pdf';txt=subprocess.check_output(['pdfinfo',str(p)],text=True);assert int(re.search(r'Pages:\s+(\d+)',txt)[1])==(2 if key=='all_cuts' else 1)
assert all(sha(out/p)==h for p,h in oldhashes.items());report['status']='PASS';report['visual_review']='pending'
(stage/'validation.json').write_text(json.dumps(report,indent=2)+'\n');shutil.copy2(__file__,stage/Path(__file__).name)
print('NUMERICAL VALIDATION PASS',flush=True)
