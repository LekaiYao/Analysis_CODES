#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""Append Prediction using the accepted three/four-curve drawing code and cached events."""
from pathlib import Path
import json, math, shutil, subprocess, tempfile, hashlib
import numpy as np
import ROOT as R
from ppref_psi2s_pthat_validation import dump
repo=Path(__file__).resolve().parents[2]
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
base=repo/'plotER/Validation'
d=dict(np.load(base/'v27_xeff48_signal_methods/artifacts/selected_data.npz'))
mc=dict(np.load(base/'v27_xeff48_signal_methods_mc/artifacts/mc_selected.npz'))
for suffix,source in [('', 'pb23_v27_optimal_splot_sideband.py'),('_mc','plot_v27_optimal_signal_methods_mc.py')]:
 target=base/('v27_xeff48_signal_methods'+suffix)
 assert not (target/'pdf/Prediction.pdf').exists()
 s=json.loads((target/'summary.json').read_text());previous=list(s['variables']);assert len(previous)==27
 protected={str(p.relative_to(target)):sha(p) for p in (target/'pdf').glob('*.pdf')}
 weights={k:d[k+'_weight'] for k in ['splot','subtraction','sideband']}
 if suffix:weights['mc']=mc['pThatreweight'].astype(float)
 names=['Prediction'];ranges={'Prediction':(s['threshold'],1.,False)};titles={'Prediction':'Prediction (XGBoost score)'}
 assert (d['Prediction']<1).all() and (mc['Prediction']<1).all()
 sourcepath=repo/'scripts/workflows'/source;text=sourcepath.read_text()
 block=text[text.index("root=R.TFile("):text.index('root.Close();dump')]+"root.Close()"
 with tempfile.TemporaryDirectory(prefix='v27_prediction_') as td:
  out=Path(td);art=out/'artifacts';art.mkdir();(out/'pdf').mkdir()
  exec(compile(block,str(sourcepath), 'exec'),globals())
  backup=target/'artifacts/before_prediction';backup.mkdir()
  for rel in ['validation_all.pdf','variable_metrics.json','summary.json','README.md','delivery_hashes.json']:
   if (target/rel).exists():shutil.copy2(target/rel,backup/rel)
  shutil.copy2(out/'pdf/Prediction.pdf',target/'pdf/Prediction.pdf')
  shutil.copy2(art/'histograms.root',target/'artifacts/prediction_histograms.root')
  metrics=json.loads((target/'variable_metrics.json').read_text());dump(target/'variable_metrics.json',records+metrics)
  s['variables']=['Prediction']+previous;s['prediction_range']=[s['threshold'],1.];s['prediction_extension_refit']=False;dump(target/'summary.json',s)
  subprocess.run(['pdfunite',*[str(target/'pdf'/f'{n}.pdf') for n in s['variables']],str(out/'validation_all.pdf')],check=True)
  shutil.copy2(out/'validation_all.pdf',target/'validation_all.pdf')
  readme=(target/'README.md').read_text().replace('27页合订本','28页合订本')
  readme+='\n新增Prediction置于第一页，展示同一选后样本的score区间('+str(s['threshold'])+',1)，15个等宽bin；原27张单图保持不变。未重新拟合。\n'
  (target/'README.md').write_text(readme)
 for rel,h in protected.items():assert sha(target/rel)==h
 dump(target/'artifacts/prediction_extension.json',dict(source_script=str(sourcepath),source_script_sha256=sha(sourcepath),original_single_pdf_hashes=protected,pdf_order=s['variables'],range=s['prediction_range'],refit=False))
 shutil.copy2(Path(__file__),target/'artifacts'/Path(__file__).name)
 print('DONE',target)
