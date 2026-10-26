from pathlib import Path
import ast,json,hashlib,shutil,subprocess,math
import numpy as np,uproot,ROOT as R
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0);R.gStyle.SetTextFont(42)
repo=Path('/eos/home-l/leyao/pbpb_work/X_analysis/Analysis_CODES');base=Path('/tmp/leyao/ppref_b_baseline_20261013');p=base/'Variables_ppRef'
source=repo/'scripts/workflows/plot_ppref_b_baseline.py'
module=ast.parse(source.read_text());fn=next(x for x in module.body if isinstance(x,ast.FunctionDef) and x.name=='plot')
labels={'Bu':'B^{+}','Bs':'B_{s}^{0}','Bd':'B^{0}'}
exec(compile(ast.Module(body=[fn],type_ignores=[]),str(source),'exec'))
pro=json.loads((p/'metadata/provenance.json').read_text());configs=json.loads((p/'metadata/variables.json').read_text())
rootfiles={'Variables_ppRef':R.TFile(str(p/'metadata/histograms.root'),'UPDATE')};changed={}
for ch,c in pro['channels'].items():
 bins_path=p/'metadata'/f'{ch}_bins.json';bins=json.loads(bins_path.read_text());bad=[n for n,x in bins.items() if sum(x['counts'])==0 and c['statistics'][n]['finite']>0];changed[ch]=bad
 if not bad:continue
 before=Path(c['source']).stat();counts={}
 for n in bad:
  cfg=configs[ch][n];s=c['statistics'][n];lo,hi=s['min'],s['max']
  if cfg['absolute_display']:lo,hi=(0. if lo<=0<=hi else min(abs(lo),abs(hi))),max(abs(lo),abs(hi))
  pad=.5 if hi==lo else (hi-lo)*.025;cfg['range']=[lo-pad,hi+pad];cfg['range_source']='full finite support: legacy range has no entries';s['underflow']=s['overflow']=0;counts[n]=np.zeros(15,dtype=np.int64)
 with uproot.open(c['source']) as f:
  for a in f[c['tree']].iterate(bad,step_size='80 MB',library='np'):
   for n,x in a.items():
    cfg=configs[ch][n];x=x[np.isfinite(x)].astype(float);x=abs(x) if cfg['absolute_display'] else x;lo,hi=cfg['range'];assert ((x>=lo)&(x<hi)).all();counts[n]+=np.histogram(x,np.linspace(lo,hi,16))[0]
 for n in bad:
  s=c['statistics'][n];assert int(sum(counts[n]))==s['finite'];bins[n]={'counts':counts[n].tolist(),'config':configs[ch][n]};plot(ch,n,configs[ch][n],counts[n],s,p/ch/'pdf'/f'{n}.pdf')
 after=Path(c['source']).stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)==(c['fingerprint']['size'],c['fingerprint']['mtime_ns'])
 bins_path.write_text(json.dumps(bins,indent=2)+'\n')
 subprocess.run(['pdfunite',*[str(p/ch/'pdf'/f'{n}.pdf') for n in c['variables']],str(p/ch/'variables_all.pdf')],check=True)
rootfiles['Variables_ppRef'].Close()
with uproot.open(p/'metadata/histograms.root') as f:
 for ch,c in pro['channels'].items():
  bins=json.loads((p/'metadata'/f'{ch}_bins.json').read_text())
  for n,x in bins.items():
   h=f[ch+'_'+n];assert np.array_equal(h.values(),x['counts']);assert np.allclose(h.variances(),x['counts'],rtol=1e-14,atol=1e-12)
pro['display_adjustments']=changed;pro['final_script_sha256']=hashlib.sha256(source.read_bytes()).hexdigest();pro['display_fix_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(p/'metadata/provenance.json').write_text(json.dumps(pro,indent=2)+'\n');(p/'metadata/variables.json').write_text(json.dumps(configs,indent=2)+'\n')
for folder in [p,base/'Bmass_ppRef']:
 shutil.copy2(source,folder/'metadata'/source.name)
shutil.copy2(__file__,p/'metadata'/Path(__file__).name)
r=p/'README.md';r.write_text(r.read_text()+'\n旧显示范围完全没有候选的变量已改为完整有限取值范围，具体列表见 provenance.json 的 display_adjustments；仅调整显示，未改变 selection。\n')
print('DISPLAY CORRECTIONS PASS',changed,flush=True)
