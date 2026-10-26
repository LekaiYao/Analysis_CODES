#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""Two-species validation, reusing the accepted ROOT 15-bin sPlot layout."""
from pathlib import Path
import argparse,json,csv,re,subprocess,shutil
import numpy as np,uproot,ROOT as R
from ppref_psi2s_pthat_validation import cdf,corr,dump,fingerprint

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);out=ap.parse_args().output.resolve();repo=Path(__file__).resolve().parents[2]
 summary=json.loads((out/'summary.json').read_text());assert summary['status']=='SPLOT_PASS';names=summary['variables'];R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0)
 def read(p,t):
  with uproot.open(p) as f:return f[t].arrays(library='np')
 d=read(out/'artifacts/data_sweights.root','ntmix_sweights');psiold={r['variable']:r for r in json.loads((repo/'plotER/validation_archive/archive/previous_validation/results/psi2s_ppref_validation/variable_metrics.json').read_text())};xold={r['variable']:r for r in csv.DictReader((repo/'plotER/validation_archive/archive/previous_validation/results/ppRef_X_r5_splot_ppref_snapshot_v1/validation/common_scalar_cdf_metrics.csv').open())}
 config={}
 for line in (repo/'plotER/validation_archive/aux.h').read_text().splitlines():
  q=re.search(r'\{"([^\"]+)\",\s*\"([^\"]+)\",\s*kNBins,\s*([^,]+),\s*([^,]+),\s*(true|false)',line)
  if q:
   try:config[q[1]]=(float(q[3]),float(q[4]),q[5]=='true',q[2].split(';')[1])
   except ValueError:pass
 R.gInterpreter.Declare('#include "'+str(repo/'plotER/validation_archive/aux.h')+'"')
 allrows={}
 for species,tree,label in [('x','ntmix_X3872','PbPb23 X(3872)'),('psi','ntmix_PSI2S','PbPb23 #psi(2S)')]:
  dest=out/species;dest.mkdir();fig=dest/'pdf';fig.mkdir();m=read(out/f'artifacts/{species}_mc.root',tree);cfg=dict(config)
  for q in R.getSignalVars(tree):cfg[str(q.expr)]=(float(q.xmin),float(q.xmax),bool(q.absVal),str(q.title).split(';')[1])
  rows=[];root=R.TFile(str(dest/'histograms.root'),'RECREATE')
  for j,n in enumerate(names):
   good=np.isfinite(d[n]);gm=np.isfinite(m[n]);x=d[n][good].astype(float);w=d[species+'_sWeight'][good];y=m[n][gm].astype(float);v=m['pThatreweight'][gm].astype(float);assert w.sum()>0 and v.sum()>0
   _,_,_,distance=cdf(x,w,y,v)
   old=psiold.get(n) if species=='psi' else xold.get(n)
   fallback=(float(old['display_min']),float(old['display_max']),False,n) if species=='psi' and old else ((float(old['support_min']),float(old['support_max']),False,n) if old else None)
   assert n in cfg or fallback is not None,n
   lo,hi,ab,title=cfg.get(n,fallback)
   if hi<=lo:lo,hi=lo-.5,hi+.5
   xx=np.abs(x) if ab else x;yy=np.abs(y) if ab else y;edges=np.linspace(lo,hi,16)
   rec={'variable':n,'cdf_distance':distance,'mc_mass_pearson':corr(y,m['Bmass'][gm],v),'data_nonfinite':int((~good).sum()),'mc_nonfinite':int((~gm).sum()),'constant_both':bool(np.ptp(x)==0 and np.ptp(y)==0),'display_min':lo,'display_max':hi,'display_abs':ab,'range_source':'aux.h' if n in cfg else 'accepted_species_baseline','bins':15};hs=[]
   for key,vals,weights in [('data',xx,w),('mc',yy,v)]:
    inside=(vals>=lo)&(vals<hi);raw=np.histogram(vals[inside],edges,weights=weights[inside])[0];err=np.sqrt(np.histogram(vals[inside],edges,weights=weights[inside]**2)[0]);norm=raw.sum()
    rec[key+'_inrange_weight_fraction']=float(norm/weights.sum());rec[key+'_hist_normalization']='inrange' if norm>0 else 'full_support_fallback'
    if norm<=0:norm=weights.sum()
    h=R.TH1D(n+'_'+key,'',15,lo,hi);h.SetDirectory(0);h.Sumw2()
    for i in range(15):h.SetBinContent(i+1,raw[i]/norm);h.SetBinError(i+1,err[i]/norm)
    hs.append(h)
   sp,mc=hs;sp.SetLineColor(R.kRed+1);sp.SetMarkerColor(R.kRed+1);sp.SetMarkerStyle(24);sp.SetLineWidth(2);mc.SetLineColor(R.kOrange+7);mc.SetLineWidth(2)
   band=mc.Clone(n+'_band');band.SetDirectory(0);band.SetFillColorAlpha(R.kOrange+7,.3);band.SetMarkerSize(0)
   c=R.TCanvas('c'+str(j),'',760,650);top=R.TPad('top'+str(j),'',0,.30,1,1);bot=R.TPad('bot'+str(j),'',0,0,1,.30)
   for pad in [top,bot]:pad.SetLeftMargin(.14);pad.SetRightMargin(.04);pad.Draw()
   top.SetBottomMargin(.01);bot.SetTopMargin(.01);bot.SetBottomMargin(.33);top.cd()
   high=max(h.GetBinContent(i)+h.GetBinError(i) for h in hs for i in range(1,16));low=min(0.,min(sp.GetBinContent(i)-sp.GetBinError(i) for i in range(1,16)))
   sp.SetMaximum(max(.1,high)*1.7);sp.SetMinimum(low*1.15);sp.GetYaxis().SetTitle('Normalized entries');sp.GetXaxis().SetLabelSize(0);sp.GetXaxis().SetTitleSize(0);sp.Draw('E');band.Draw('E2 SAME');mc.Draw('HIST SAME');sp.Draw('E SAME')
   leg=R.TLegend(.62,.77,.94,.90);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.04);leg.AddEntry(mc,'Prompt MC','l');leg.AddEntry(sp,'sPlot DATA','lep');leg.Draw()
   lab=R.TLatex();lab.SetNDC();lab.SetTextSize(.052);lab.DrawLatex(.18,.87,label);lab.SetTextSize(.034);lab.DrawLatex(.18,.80,'v27 v1, X efficiency 50%');lab.DrawLatex(.18,.735,'CDF D = %.3f'%distance);lab.SetTextSize(.029);lab.DrawLatex(.18,.67,'MC: pThatreweight; 15 bins')
   if min(rec['data_inrange_weight_fraction'],rec['mc_inrange_weight_fraction'])<.95:lab.DrawLatex(.18,.61,'In-range fraction DATA/MC: %.3f / %.3f'%(rec['data_inrange_weight_fraction'],rec['mc_inrange_weight_fraction']))
   bot.cd();axis=R.TH1D(n+'_ratio_axis','',15,lo,hi);axis.SetDirectory(0);axis.SetMinimum(0);axis.SetMaximum(5);axis.GetYaxis().SetTitle('sPlot / MC');axis.GetYaxis().SetTitleSize(.09);axis.GetYaxis().SetLabelSize(.08);axis.GetYaxis().SetTitleOffset(.7);axis.GetYaxis().SetNdivisions(304);axis.GetXaxis().SetTitle(title);axis.GetXaxis().SetTitleSize(.11);axis.GetXaxis().SetLabelSize(.10);axis.Draw('AXIS')
   g=R.TGraphErrors();g.SetName(n+'_ratio');clip=0
   for i in range(1,16):
    den=mc.GetBinContent(i);num=sp.GetBinContent(i)
    if den<=0:continue
    val=num/den;err=np.hypot(sp.GetBinError(i)/den,num*mc.GetBinError(i)/den**2);k=g.GetN();g.SetPoint(k,mc.GetBinCenter(i),val);g.SetPointError(k,0,float(err));clip+=int(val<0 or val>5)
   g.SetMarkerStyle(20);g.SetMarkerColor(R.kOrange+7);g.SetLineColor(R.kOrange+7);g.Draw('P SAME');line=R.TLine(lo,1,hi,1);line.SetLineStyle(2);line.Draw();rec['ratio_points_outside_0_5']=clip
   if clip:lab.SetTextSize(.075);lab.DrawLatex(.18,.82,f'{clip} ratio points outside [0,5]')
   c.SaveAs(str(fig/f'{n}.pdf'));root.cd()
   for obj in [sp,mc,g]:obj.Write()
   c.Close();rows.append(rec)
  root.Close();subprocess.run(['pdfunite',*[str(fig/f'{n}.pdf') for n in names],str(out/('x_validation_all.pdf' if species=='x' else 'psi2s_validation_all.pdf'))],check=True);dump(dest/'variable_metrics.json',rows)
  with (dest/'variable_metrics.csv').open('w') as f:
   writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
  allrows[species]=rows;print('PLOTTED',species,len(rows),flush=True)
 shutil.copy2(__file__,out/'artifacts'/Path(__file__).name);shutil.copy2(repo/'scripts/workflows/ppref_psi2s_prompt_plots.py',out/'artifacts/accepted_plot_template.py')
 dump(out/'plot_metadata.json',{'template':fingerprint(repo/'scripts/workflows/ppref_psi2s_prompt_plots.py'),'script':fingerprint(Path(__file__)),'aux':fingerprint(repo/'plotER/validation_archive/aux.h'),'fallback_psi':fingerprint(repo/'plotER/validation_archive/archive/previous_validation/results/psi2s_ppref_validation/variable_metrics.json'),'fallback_x':fingerprint(repo/'plotER/validation_archive/archive/previous_validation/results/ppRef_X_r5_splot_ppref_snapshot_v1/validation/common_scalar_cdf_metrics.csv'),'variables_per_species':len(names),'bins':15,'ratio_display':[0,5],'negative_sweights_retained':True,'cdf_full_finite_raw_support':True})
 summary['status']='PLOTS_COMPLETE_PENDING_AUDIT';dump(out/'summary.json',summary)
if __name__=='__main__':main()
