#!/usr/bin/env python3
"""Render both FOM coordinates from frozen scan results; no fits or point selection."""
from pathlib import Path
import json,hashlib,argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot(out):
 out=Path(out);rows=json.loads((out/'scan.json').read_text());selected=json.loads((out/'selected_point.json').read_text())
 assert hashlib.sha256((out/'scan.json').read_bytes()).hexdigest()==selected['scan_sha256']
 best=next(r for r in rows if r['xeff']==selected['xeff']);assert best['score_threshold']==selected['score_threshold']
 for key,name,label in [('xeff','fom_scan.pdf','Weighted MC efficiency [%]'),('score_threshold','fom_score.pdf','Score threshold (Prediction > threshold)')]:
  ordered=sorted(rows,key=lambda r:r[key]);x=[r[key] for r in ordered]
  fig,ax=plt.subplots(figsize=(7,5))
  ax.errorbar(x,[r['fom'] for r in ordered],yerr=[r['fom_background_error'] for r in ordered],fmt='o-',ms=3,capsize=2,label='Reference yield 78.79; background error only')
  ax.plot(x,[r['fom_low_reference'] for r in ordered],'--',label='Reference yield -1 sigma')
  ax.plot(x,[r['fom_high_reference'] for r in ordered],'--',label='Reference yield +1 sigma')
  ax.axvline(best[key],color='gray',ls=':',lw=1)
  ax.plot(best[key],best['fom'],'*',color='crimson',ms=12,zorder=5)
  ax.annotate('Maximum: xeff = 36%\nscore > 0.850995',xy=(best[key],best['fom']),xytext=(.5,.68),textcoords='axes fraction',ha='center',fontsize=9,arrowprops={'arrowstyle':'-','color':'crimson'})
  ax.set(xlabel=label,ylabel=r'$S/\sqrt{S+B}$',title='PbPb23 v27: Run-2 counting FOM')
  ax.legend(fontsize=8,loc='center',bbox_to_anchor=(.5,.34),framealpha=1.)
  fig.tight_layout();fig.savefig(out/name);plt.close(fig)
 meta={'scan_sha256':selected['scan_sha256'],'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'no_refit':True,'figures':{n:hashlib.sha256((out/n).read_bytes()).hexdigest() for n in ['fom_scan.pdf','fom_score.pdf']},'error_definition':'background statistical error only; reference +/-1 sigma are sensitivity curves, not confidence bands'}
 (out/'plot_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
if __name__=='__main__':
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('output',type=Path);plot(ap.parse_args().output)
