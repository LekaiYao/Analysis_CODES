#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""Replot frozen Punzi rows without refitting or changing grid/selection."""
from pathlib import Path
import json,hashlib,argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def plot(out):
 out=Path(out)
 contract=json.loads((out/'contract.json').read_text());label_tag=contract.get('plot_label','PbPb23 v27 v1')
 for grid in ['xeff','score']:
  rows=json.loads((out/f'{grid}_scan.json').read_text())
  valid=[r for r in rows if r['status']=='PASS']
  for coord,label in [('xeff_actual_percent','Weighted MC efficiency [%]'),('score_threshold','Score threshold (Prediction > threshold)')]:
   ordered=sorted(valid,key=lambda r:r[coord]);fig,ax=plt.subplots(figsize=(7,5))
   for a in [3,5]:
    key=f'punzi_a{a}';line=ax.errorbar([r[coord] for r in ordered],[r[key] for r in ordered],yerr=[r[key+'_background_error'] for r in ordered],fmt='o-',ms=3,capsize=2,label=f'a = {a}; background error only')
    best=max(valid,key=lambda r:r[key]);ax.plot(best[coord],best[key],'*',ms=12,color=line[0].get_color());ax.axvline(best[coord],ls=':',lw=1,color=line[0].get_color())
   best=max(valid,key=lambda r:r['punzi_a3'])
   other=max(valid,key=lambda r:r['punzi_a5'])
   ax.text(.5,.18,f"a=3 max: xeff={best['xeff_actual_percent']:.3f}%, score>{best['score_threshold']:.6f}\na=5 max: xeff={other['xeff_actual_percent']:.3f}%, score>{other['score_threshold']:.6f}",transform=ax.transAxes,ha='center',va='top',fontsize=9)
   ax.set(xlabel=label,ylabel='Punzi FOM',title=f'{label_tag}: {grid} grid');ax.legend(fontsize=8,loc='lower center')
   fig.tight_layout();fig.savefig(out/f'{grid}_scan_vs_{"xeff" if coord=="xeff_actual_percent" else "score"}.pdf');plt.close(fig)
 meta={'scan_hashes':{g:hashlib.sha256((out/f'{g}_scan.json').read_bytes()).hexdigest() for g in ['xeff','score']},'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'no_refit':True,'pdf_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob('*.pdf')}}
 (out/'plot_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('output',type=Path);plot(ap.parse_args().output)
