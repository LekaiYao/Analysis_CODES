#!/usr/bin/env python3
"""ppRef three-channel mass spectra: fiducial + 2D4; ML sidebands drawn only."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,shutil,subprocess,re,copy
import numpy as np
import uproot
import ROOT as R
from ppref_b_plot_style import plot
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0);R.gStyle.SetTextFont(42)
assert R.gROOT.GetVersion()=='6.32.02'
repo=Path(__file__).resolve().parents[2];out=repo/'plotER/results/Bmass_ppRef';stage=Path('/tmp/leyao/ppref_b_mass_fid_2d4');assert not stage.exists();stage.mkdir()
key='Bpt7p5_By2p4_2D4';sel='Bpt > 7.5 && abs(By) < 2.4 && Bnorm_svpvDistance_2D > 4'
pro=json.loads((out/'metadata/provenance.json').read_text());cfgs=json.loads((out/'metadata/variables.json').read_text())
sidebands={'Bu':[[5.,5.2],[5.36,5.56]],'Bd':[[5.,5.2],[5.36,5.56]],'Bs':[[5.10,5.29],[5.45,5.64]]}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
oldhash={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()};(stage/'before_hashes.json').write_text(json.dumps(oldhash,indent=2)+'\n')
root=R.TFile(str(stage/'histograms.root'),'RECREATE');report={'created_UTC':datetime.now(timezone.utc).isoformat(),'selection':sel,'name':key,'display_expression':'Bmass > 5 && Bmass < 6','sideband_is_selection':False,'sideband_origin':'Bu ppRef training; Bd/Bs channel definitions from PbPb training, shown on ppRef by user request','sidebands':sidebands,'weights':'unit candidate counts','source_SHA_recomputed':False,'channels':{},'script_sha256':sha(Path(__file__)),'style_sha256':sha(Path(__file__).with_name('ppref_b_plot_style.py'))}
for ch in ['Bu','Bs','Bd']:
 c=pro['channels'][ch];src=Path(c['source']);before=src.stat();assert (before.st_size,before.st_mtime_ns)==(c['fingerprint']['size'],c['fingerprint']['mtime_ns'])
 tag=f'{ch}_{"pp24" if ch=="Bu" else "pb23"}_v5_fid3_3v1_rw0_xgb_v1';mp=repo.parent/'XGBoost/output/models'/tag/'run_metadata.json';md=json.loads(mp.read_text());match=re.search(r'\(\(\(Bmass > ([\d.]+)\) and \(Bmass < ([\d.]+)\)\) or \(\(Bmass > ([\d.]+)\) and \(Bmass < ([\d.]+)\)\)\)',md['background_selection']);assert match,md['background_selection'];assert list(map(float,match.groups()))==sum(sidebands[ch],[])
 s={'total':0,'finite':0,'nonfinite':0,'min':None,'max':None,'underflow':0,'overflow':0};counts=np.zeros(100,dtype=np.int64);total=0
 with uproot.open(src) as f:
  for a in f[c['tree']].iterate(['Bmass','Bpt','By','Bnorm_svpvDistance_2D'],step_size='80 MB',library='np'):
   pt=a['Bpt'].astype(float);y=a['By'].astype(float);d=a['Bnorm_svpvDistance_2D'].astype(float);assert np.isfinite(pt).all() and np.isfinite(y).all() and np.isfinite(d).all();mask=(pt>7.5)&(abs(y)<2.4)&(d>4);x=a['Bmass'][mask].astype(float);total+=len(pt);s['total']+=len(x);ok=np.isfinite(x);s['nonfinite']+=int((~ok).sum());x=x[ok];s['finite']+=len(x)
   if len(x):s['min']=float(x.min()) if s['min'] is None else min(s['min'],float(x.min()));s['max']=float(x.max()) if s['max'] is None else max(s['max'],float(x.max()))
   counts+=np.histogram(x[(x>5)&(x<6)],np.linspace(5,6,101))[0];s['underflow']+=int((x<=5).sum());s['overflow']+=int((x>=6).sum())
 assert total==c['entries'];assert s['total']=={'Bu':273455,'Bs':12039,'Bd':89570}[ch];assert sum(counts)+s['underflow']+s['overflow']==s['finite']
 df=R.RDataFrame(c['tree'],str(src)).Filter(sel);n=df.Count();hp=df.Filter('Bmass > 5 && Bmass < 6').Histo1D((ch+'_audit','',100,5.,6.),'Bmass');low=df.Filter('Bmass <= 5').Count();high=df.Filter('Bmass >= 6').Count();h=hp.GetValue();assert n.GetValue()==s['total'];assert np.array_equal(counts,[h.GetBinContent(i) for i in range(1,101)]);assert np.allclose(counts,[h.GetBinError(i)**2 for i in range(1,101)]);assert low.GetValue()==s['underflow'] and high.GetValue()==s['overflow']
 cfg=copy.deepcopy(cfgs[ch]['Bmass']);assert cfg['range']==[5,6] and cfg['bins']==100;cfg.update(sidebands=sidebands[ch],fiducial_label='p_{T} > 7.5 GeV/c, |y| < 2.4')
 pdf=stage/f'{ch}_{key}.pdf';plot(ch,'Bmass_'+key,cfg,counts,s,pdf,{'Bmass_ppRef':root},mass=True,selection_label='Bnorm_svpvDistance_2D > 4')
 subprocess.run(['pdfunite',str(out/ch/'pdf/baseline.pdf'),str(out/ch/'pdf/2D4.pdf'),str(pdf),str(stage/f'{ch}_all_cuts.pdf')],check=True)
 after=src.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
 report['channels'][ch]={'source':str(src),'tree':c['tree'],'fingerprint':c['fingerprint'],'input_entries':total,'statistics':s,'config':cfg,'bin_counts':counts.tolist(),'in_range':int(sum(counts)),'ML_metadata':str(mp),'ML_metadata_sha256':sha(mp),'ROOT_numpy_selection_bins_errors':'PASS'}
 print(ch,s['total'],'selected;',int(sum(counts)),'in range PASS',flush=True)
root.Close()
with uproot.open(stage/'histograms.root') as f:
 for ch,c in report['channels'].items():h=f[ch+'_Bmass_'+key];assert np.array_equal(h.values(),c['bin_counts']);assert np.allclose(h.variances(),c['bin_counts'])
for ch in report['channels']:
 for suffix,pages in [(key,1),('all_cuts',3)]:
  pdf=stage/f'{ch}_{suffix}.pdf';info=subprocess.check_output(['pdfinfo',str(pdf)],text=True);assert int(re.search(r'Pages:\s+(\d+)',info)[1])==pages
assert all(sha(out/p)==h for p,h in oldhash.items());report['status']='PASS';report['visual_review']='pending'
(stage/'provenance.json').write_text(json.dumps(report,indent=2)+'\n')
for p in [Path(__file__),Path(__file__).with_name('ppref_b_plot_style.py')]:shutil.copy2(p,stage/p.name)
print('NUMERICAL VALIDATION PASS',flush=True)
