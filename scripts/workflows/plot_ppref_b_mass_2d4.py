#!/usr/bin/env python3
"""ppRef raw DATA with the sole extra cut Bnorm_svpvDistance_2D > 4."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib,json,shutil,subprocess,re
import numpy as np
import uproot
import ROOT as R
from ppref_b_plot_style import plot
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0);R.gStyle.SetTextFont(42)
assert R.gROOT.GetVersion()=='6.32.02'
repo=Path(__file__).resolve().parents[2];out=repo/'plotER/results/Bmass_ppRef'
stage=Path('/tmp/leyao/ppref_b_mass_2d4');assert not stage.exists();stage.mkdir()
pro=json.loads((out/'metadata/provenance.json').read_text());cfgs=json.loads((out/'metadata/variables.json').read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
protected={str(p):sha(p) for p in out.rglob('*') if p.is_file()}
rootfiles={'Bmass_ppRef':R.TFile(str(stage/'histograms.root'),'RECREATE')}
report={'selection':'Bnorm_svpvDistance_2D > 4','name':'2D4','created_UTC':datetime.now(timezone.utc).isoformat(),'ROOT_version':R.gROOT.GetVersion(),'additional_cuts':['Bnorm_svpvDistance_2D > 4'],'weights':'unit candidate counts','channels':{},'script_sha256':sha(Path(__file__)),'style_sha256':sha(Path(__file__).with_name('ppref_b_plot_style.py'))}
for ch in ['Bu','Bs','Bd']:
 c=pro['channels'][ch];src=Path(c['source']);before=src.stat();assert (before.st_size,before.st_mtime_ns)==(c['fingerprint']['size'],c['fingerprint']['mtime_ns'])
 cfg=cfgs[ch]['Bmass'];lo,hi=cfg['range'];bins=cfg['bins'];counts=np.zeros(bins,dtype=np.int64)
 s={'total':0,'finite':0,'nonfinite':0,'min':None,'max':None,'underflow':0,'overflow':0};total=0;bad_cut=0
 with uproot.open(src) as f:
  for a in f[c['tree']].iterate(['Bmass','Bnorm_svpvDistance_2D'],step_size='80 MB',library='np'):
   cut=a['Bnorm_svpvDistance_2D'];assert np.isfinite(cut).all();bad_cut+=int((~np.isfinite(cut)).sum());mask=cut>4;total+=len(cut)
   mass=a['Bmass'][mask].astype(float);s['total']+=len(mass);ok=np.isfinite(mass);s['nonfinite']+=int((~ok).sum());mass=mass[ok];s['finite']+=len(mass)
   if len(mass):s['min']=float(mass.min()) if s['min'] is None else min(s['min'],float(mass.min()));s['max']=float(mass.max()) if s['max'] is None else max(s['max'],float(mass.max()))
   inside=(mass>lo)&(mass<hi);counts+=np.histogram(mass[inside],np.linspace(lo,hi,bins+1))[0];s['underflow']+=int((mass<=lo).sum());s['overflow']+=int((mass>=hi).sum())
 assert total==c['entries'];assert sum(counts)+s['underflow']+s['overflow']==s['finite']
 df=R.RDataFrame(c['tree'],str(src)).Filter(report['selection']);n=df.Count();hp=df.Filter('Bmass > 5 && Bmass < 6').Histo1D((ch+'_audit','',bins,lo,hi),'Bmass');h=hp.GetValue();assert n.GetValue()==s['total']
 assert np.array_equal(counts,[h.GetBinContent(i) for i in range(1,bins+1)])
 assert df.Filter('Bmass <= 5').Count().GetValue()==s['underflow'] and df.Filter('Bmass >= 6').Count().GetValue()==s['overflow']
 assert np.allclose([h.GetBinError(i)**2 for i in range(1,bins+1)],counts)
 pdf=stage/(ch+'_2D4.pdf');plot(ch,'Bmass_2D4',cfg,counts,s,pdf,rootfiles,mass=True,selection_label='Bnorm_svpvDistance_2D > 4')
 subprocess.run(['pdfunite',str(out/ch/'pdf/baseline.pdf'),str(pdf),str(stage/(ch+'_all_cuts.pdf'))],check=True)
 after=src.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
 report['channels'][ch]={'source':str(src),'tree':c['tree'],'fingerprint':c['fingerprint'],'input_entries':total,'statistics':s,'config':cfg,'bin_counts':counts.tolist(),'ROOT_numpy_selection_bins_errors':'PASS'}
 print(ch,s['total'],'selected;',int(sum(counts)),'in range',flush=True)
rootfiles['Bmass_ppRef'].Close()
with uproot.open(stage/'histograms.root') as f:
 for ch,c in report['channels'].items():
  h=f[ch+'_Bmass_2D4'];assert np.array_equal(h.values(),c['bin_counts']);assert np.allclose(h.variances(),c['bin_counts'])
for ch in report['channels']:
 for suffix,pages in [('2D4',1),('all_cuts',2)]:
  pdf=stage/f'{ch}_{suffix}.pdf';info=subprocess.check_output(['pdfinfo',str(pdf)],text=True);assert int(re.search(r'Pages:\s+(\d+)',info)[1])==pages
 assert 'Bnorm_svpvDistance_2D > 4' in subprocess.check_output(['pdftotext',str(stage/f'{ch}_2D4.pdf'),'-'],text=True)
assert all(sha(Path(p))==h for p,h in protected.items())
report['status']='PASS';report['baseline_and_existing_artifacts_unchanged']=True;report['source_SHA_recomputed']=False
(stage/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')
for p in [Path(__file__),Path(__file__).with_name('ppref_b_plot_style.py')]:shutil.copy2(p,stage/p.name)
print('NUMERICAL VALIDATION PASS',flush=True)
