#!/usr/bin/env python3
"""PbPb23 B DATA mass: fiducial, fiducial+2D4; grey ML-sideband guides."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,shutil,subprocess,re,copy,time
import numpy as np
import uproot
import ROOT as R
from ppref_b_plot_style import plot
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0);R.gStyle.SetTextFont(42)
assert R.gROOT.GetVersion()=='6.32.02'
repo=Path(__file__).resolve().parents[2];out=Path('/tmp/leyao/Bmass_PbPb23');assert not out.exists();(out/'metadata').mkdir(parents=True)
invpath=repo.parent/'XGBoost/output/internal/code_validation/b_samples_request_20261013/inventory.json';inventory=json.loads(invpath.read_text())
basecfg=json.loads((repo/'plotER/configs/b_inputs/ppref_mass.json').read_text())
selections={'Bpt7p5_By2p4':'Bpt > 7.5 && abs(By) < 2.4','Bpt7p5_By2p4_2D4':'Bpt > 7.5 && abs(By) < 2.4 && Bnorm_svpvDistance_2D > 4'}
sb={'Bu':[[5.,5.2],[5.36,5.56]],'Bd':[[5.,5.2],[5.36,5.56]],'Bs':[[5.1,5.29],[5.45,5.64]]}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
report={'status':'running','created_UTC':datetime.now(timezone.utc).isoformat(),'system':'PbPb23','selections':selections,'display_expression':'Bmass > 5 && Bmass < 6','bins':100,'bin_width_GeV':.01,'sidebands':sb,'sideband_color':'grey (ROOT kGray+2)','sideband_is_selection':False,'weights':'unit candidate counts','ROOT_version':R.gROOT.GetVersion(),'source_SHA_recomputed':False,'numpy_io':'XRootD direct to the same EOS path; local metadata checks','io_validation':'Bd first 1000 candidates in all four columns exactly matched local read; first million direct read 2.14 seconds','inventory_sha256':sha(invpath),'script_sha256':sha(Path(__file__)),'style_sha256':sha(Path(__file__).with_name('ppref_b_plot_style.py')),'channels':{}}
root=R.TFile(str(out/'metadata/histograms.root'),'RECREATE')
for ch in ['Bu','Bs','Bd']:
 record=next(x for x in inventory['inputs'] if x['channel']==ch and x['system']=='PbPb23' and x['kind']=='DATA');src=Path(record['path']);before=src.stat();assert (before.st_size,before.st_mtime_ns)==(record['size'],record['mtime_ns'])
 mp=repo.parent/f'XGBoost/output/models/{ch}_pb23_v5_fid3_3v1_rw0_xgb_v1/run_metadata.json';md=json.loads(mp.read_text());m=re.search(r'\(\(\(Bmass > ([\d.]+)\) and \(Bmass < ([\d.]+)\)\) or \(\(Bmass > ([\d.]+)\) and \(Bmass < ([\d.]+)\)\)\)',md['background_selection']);assert m and list(map(float,m.groups()))==sum(sb[ch],[])
 stats={k:dict(total=0,finite=0,nonfinite=0,min=None,max=None,underflow=0,overflow=0,equal_lower_edge=0,equal_upper_edge=0) for k in selections};hist={k:np.zeros(100,dtype=np.int64) for k in selections};nonfinite={n:0 for n in ['Bpt','By','Bnorm_svpvDistance_2D']};total=0;last=time.monotonic()
 print(ch,'reading',record['entries'],'candidates',flush=True)
 with uproot.open('root://eosuser.cern.ch/'+str(src)) as f:
  t=f[record['tree']];assert t.num_entries==record['entries']
  for a in t.iterate(['Bmass','Bpt','By','Bnorm_svpvDistance_2D'],step_size='100 MB',library='np'):
   pt=a['Bpt'].astype(float);y=a['By'].astype(float);d=a['Bnorm_svpvDistance_2D'].astype(float);base=(pt>7.5)&(abs(y)<2.4);total+=len(pt)
   for n in nonfinite:nonfinite[n]+=int((~np.isfinite(a[n])).sum())
   for key,mask in zip(selections,[base,base&(d>4)]):
    x=a['Bmass'][mask].astype(float);s=stats[key];s['total']+=len(x);ok=np.isfinite(x);s['nonfinite']+=int((~ok).sum());x=x[ok];s['finite']+=len(x)
    if len(x):s['min']=float(x.min()) if s['min'] is None else min(s['min'],float(x.min()));s['max']=float(x.max()) if s['max'] is None else max(s['max'],float(x.max()))
    hist[key]+=np.histogram(x[(x>5)&(x<6)],np.linspace(5,6,101))[0];s['underflow']+=int((x<=5).sum());s['overflow']+=int((x>=6).sum());s['equal_lower_edge']+=int((x==5).sum());s['equal_upper_edge']+=int((x==6).sum())
   if time.monotonic()-last>25:print(ch,'read',total,'/',record['entries'],flush=True);last=time.monotonic()
 assert total==record['entries'];dump(out/'metadata'/f'{ch}_numpy_checkpoint.json',{'source':str(src),'fingerprint':{'size':before.st_size,'mtime_ns':before.st_mtime_ns},'selections':selections,'stats':stats,'histograms':{k:v.tolist() for k,v in hist.items()},'selection_nonfinite':nonfinite})
 print(ch,'NumPy complete; independent ROOT pass',flush=True)
 df=R.RDataFrame(record['tree'],str(src));actions={}
 for key,sel in selections.items():
  filtered=df.Filter(sel);actions[key]={'count':filtered.Count(),'low':filtered.Filter('Bmass <= 5').Count(),'high':filtered.Filter('Bmass >= 6').Count(),'hist':filtered.Filter('Bmass > 5 && Bmass < 6').Histo1D((ch+'_'+key+'_audit','',100,5.,6.),'Bmass')}
 report['channels'][ch]={'source':str(src),'tree':record['tree'],'fingerprint':{'size':before.st_size,'mtime_ns':before.st_mtime_ns},'input_candidates':total,'selection_nonfinite':nonfinite,'ML_metadata':str(mp),'ML_metadata_sha256':sha(mp),'results':{}}
 (out/ch/'pdf').mkdir(parents=True)
 for key,ac in actions.items():
  s=stats[key];h=ac['hist'].GetValue();assert ac['count'].GetValue()==s['total'];assert np.array_equal(hist[key],[h.GetBinContent(i) for i in range(1,101)]);assert np.allclose(hist[key],[h.GetBinError(i)**2 for i in range(1,101)]);assert ac['low'].GetValue()==s['underflow'] and ac['high'].GetValue()==s['overflow'];assert sum(hist[key])+s['underflow']+s['overflow']==s['finite']
  cfg=copy.deepcopy(basecfg[ch]['Bmass']);assert cfg['range']==[5,6] and cfg['bins']==100;cfg.update(system_label='PbPb 2023',sidebands=sb[ch],sideband_color=int(R.kGray+2),fiducial_label='p_{T} > 7.5 GeV/c, |y| < 2.4' if key.endswith('2D4') else '')
  plot(ch,'Bmass_'+key,cfg,hist[key],s,out/ch/'pdf'/f'{key}.pdf',{'Bmass_ppRef':root},mass=True,selection_label='Bnorm_svpvDistance_2D > 4' if key.endswith('2D4') else 'p_{T} > 7.5 GeV/c, |y| < 2.4')
  report['channels'][ch]['results'][key]={'statistics':s,'bin_counts':hist[key].tolist(),'in_range':int(sum(hist[key])),'config':cfg,'ROOT_numpy_selection_bins_errors':'PASS'}
  print(ch,key,s['total'],'selected;',int(sum(hist[key])),'in range PASS',flush=True)
 subprocess.run(['pdfunite',*[str(out/ch/'pdf'/f'{key}.pdf') for key in selections],str(out/ch/'all_cuts.pdf')],check=True)
 after=src.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns);dump(out/'metadata/provenance.json',report)
root.Close()
with uproot.open(out/'metadata/histograms.root') as f:
 for ch,c in report['channels'].items():
  for key,d in c['results'].items():h=f[ch+'_Bmass_'+key];assert np.array_equal(h.values(),d['bin_counts']);assert np.allclose(h.variances(),d['bin_counts'])
for ch in report['channels']:
 for pdf in list((out/ch/'pdf').glob('*.pdf'))+[out/ch/'all_cuts.pdf']:
  info=subprocess.check_output(['pdfinfo',str(pdf)],text=True);assert int(re.search(r'Pages:\s+(\d+)',info)[1])==(2 if pdf.name=='all_cuts.pdf' else 1)
  text=subprocess.check_output(['pdftotext',str(pdf),'-'],text=True);assert 'PbPb 2023' in text and 'ML sideband boundaries' in text and 'not applied as a cut' in text
report['status']='PASS';report['visual_review']='pending';dump(out/'metadata/provenance.json',report);dump(out/'metadata/cuts.json',selections)
for p in [Path(__file__),Path(__file__).with_name('ppref_b_plot_style.py')]:shutil.copy2(p,out/'metadata'/p.name)
s='# Bmass_PbPb23\n\n原始PbPb23 DATA，分别施加基础fiducial（Bpt > 7.5 GeV、abs(By) < 2.4），以及再加Bnorm_svpvDistance_2D > 4。未加其它channel pre-cut、中心度cut或ML score cut。\n\n严格5 < Bmass < 6 GeV显示，100个10 MeV bins；N为全质量范围选后候选数，in range为图内候选数。误差sqrt(N)，无拟合或减背景。\n\n灰色虚线为ML侧带边界，仅示意、不筛选：Bu/Bd为(5.00,5.20)+(5.36,5.56)，Bs为(5.10,5.29)+(5.45,5.64)，单位GeV。边界已与各道pb23训练metadata核对。\n\n| Channel | 基础fiducial | 再加2D4 | 合订本 |\n|---|---|---|---|\n'
for ch in ['Bu','Bs','Bd']:s+=f'| {ch} | [Bpt7p5_By2p4]({ch}/pdf/Bpt7p5_By2p4.pdf) | [Bpt7p5_By2p4_2D4]({ch}/pdf/Bpt7p5_By2p4_2D4.pdf) | [all_cuts]({ch}/all_cuts.pdf) |\n'
s+='\n[输入、计数、分箱与验证](metadata/provenance.json)。源文件size/mtime与inventory及读前读后一致，未重算大文件SHA。仅排版重绘可复用histograms.root；selection、输入、分箱改变需重新求值。\n';(out/'README.md').write_text(s)
print('ALL SIX PLOTS NUMERICALLY VALIDATED',flush=True)
