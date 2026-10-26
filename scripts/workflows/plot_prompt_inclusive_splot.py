#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""Accepted normalized 15-bin layout; account for shared events in ratio errors."""
from pathlib import Path
import argparse,json,subprocess,shutil,re
import numpy as np
import uproot,ROOT as R
ap=argparse.ArgumentParser();ap.add_argument('--prompt',required=True,type=Path);a=ap.parse_args();out=a.prompt.resolve();repo=Path(__file__).resolve().parents[2];inc=repo/'plotER/validation_archive/results/v27_xeff50_splot'
assert R.gROOT.GetVersion()=='6.32.02';R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0)
def read(p):
 with uproot.open(p) as f:return f['ntmix_sweights'].arrays(library='np')
pd=read(out/'artifacts/data_sweights.root');di=read(inc/'artifacts/data_sweights.root');idx=np.searchsorted(di['source_entry'],pd['source_entry']);assert np.array_equal(di['source_entry'][idx],pd['source_entry'])
assert np.array_equal(di['source_entry'][di['BLxy']<.01],pd['source_entry'])
config={}
for line in (repo/'plotER/validation_archive/aux.h').read_text().splitlines():
 m=re.search(r'\{"([^"\n]+)",\s*"([^"\n]+)",\s*kNBins,',line)
 if m:config[m[1]]=m[2].split(';')[1]
R.gInterpreter.Declare('#include "'+str(repo/'plotER/validation_archive/aux.h')+'"')
report={}
for species in ['x','psi']:
 dest=out/species/'prompt_vs_inclusive';dest.mkdir();fig=dest/'pdf';fig.mkdir()
 pf=R.TFile.Open(str(out/species/'histograms.root'));fi=R.TFile.Open(str(inc/species/'histograms.root'));rf=R.TFile(str(dest/'histograms.root'),'RECREATE')
 metrics=json.loads((out/species/'variable_metrics.json').read_text());names=[r['variable'] for r in metrics];rows=[];cfg=dict(config)
 for v in R.getSignalVars('ntmix_X3872' if species=='x' else 'ntmix_PSI2S'):cfg[str(v.expr)]=str(v.title).split(';')[1]
 for j,row in enumerate(metrics):
  n=row['variable'];ph=pf.Get(n+'_data').Clone(n+'_prompt');ih=fi.Get(n+'_data').Clone(n+'_inclusive');ph.SetDirectory(0);ih.SetDirectory(0)
  lo,hi=row['display_min'],row['display_max'];edges=np.linspace(lo,hi,16);assert ih.GetXaxis().GetXmin()==lo and ih.GetXaxis().GetXmax()==hi
  x=np.abs(pd[n]) if row['display_abs'] else pd[n];y=np.abs(di[n]) if row['display_abs'] else di[n]
  good=np.isfinite(x);inside=good&(x>=lo)&(x<hi);gi=np.isfinite(y);ii=gi&(y>=lo)&(y<hi)
  wp=pd[species+'_sWeight'];wi=di[species+'_sWeight'];np_=wp[inside].sum();ni=wi[ii].sum()
  if np_<=0:np_=wp[good].sum()
  if ni<=0:ni=wi[gi].sum()
  cross=np.histogram(x[inside],edges,weights=(wp*wi[idx])[inside])[0]/(np_*ni)
  for h,col,marker in [(ph,R.kRed+1,24),(ih,R.kBlue+1,20)]:h.SetLineColor(col);h.SetMarkerColor(col);h.SetMarkerStyle(marker);h.SetLineWidth(2)
  c=R.TCanvas('c'+str(j),'',760,650);top=R.TPad('t'+str(j),'',0,.30,1,1);bot=R.TPad('b'+str(j),'',0,0,1,.30)
  for pad in [top,bot]:pad.SetLeftMargin(.14);pad.SetRightMargin(.04);pad.Draw()
  top.SetBottomMargin(.01);bot.SetTopMargin(.01);bot.SetBottomMargin(.33);top.cd()
  high=max(h.GetBinContent(i)+h.GetBinError(i) for h in [ph,ih] for i in range(1,16));low=min(0.,min(h.GetBinContent(i)-h.GetBinError(i) for h in [ph,ih] for i in range(1,16)))
  ph.SetMaximum(max(.1,high)*1.7);ph.SetMinimum(low*1.15);ph.GetYaxis().SetTitle('Normalized entries');ph.GetXaxis().SetLabelSize(0);ph.Draw('E');ih.Draw('E SAME');ph.Draw('E SAME')
  leg=R.TLegend(.58,.77,.94,.90);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.034);leg.AddEntry(ph,'Prompt-enriched sPlot','lep');leg.AddEntry(ih,'Inclusive sPlot','lep');leg.Draw()
  lab=R.TLatex();lab.SetNDC();lab.SetTextSize(.047);lab.DrawLatex(.18,.87,'PbPb23 '+('X(3872)' if species=='x' else '#psi(2S)'));lab.SetTextSize(.031);lab.DrawLatex(.18,.80,'v27 v1, X efficiency 50%');lab.DrawLatex(.18,.735,'Prompt-enriched: L_{xy} < 0.01 cm');lab.SetTextSize(.028);lab.DrawLatex(.18,.67,'Shared-event covariance in ratio errors')
  bot.cd();axis=R.TH1D(n+'_axis','',15,lo,hi);axis.SetDirectory(0);axis.SetMinimum(0);axis.SetMaximum(5);axis.GetYaxis().SetTitle('Prompt / Incl.');axis.GetYaxis().SetTitleSize(.09);axis.GetYaxis().SetLabelSize(.08);axis.GetYaxis().SetTitleOffset(.7);axis.GetYaxis().SetNdivisions(304);axis.GetXaxis().SetTitle(cfg.get(n,n));axis.GetXaxis().SetTitleSize(.11);axis.GetXaxis().SetLabelSize(.10);axis.Draw('AXIS')
  g=R.TGraphErrors();g.SetName(n+'_ratio');clips=0;ratios=[]
  for i in range(1,16):
   den=ih.GetBinContent(i);num=ph.GetBinContent(i)
   if den<=0:continue
   val=num/den;variance=(ph.GetBinError(i)/den)**2+(num*ih.GetBinError(i)/den**2)**2-2*num*cross[i-1]/den**3
   assert variance>=-1e-9,(n,variance)
   err=np.sqrt(max(0.,variance));k=g.GetN();g.SetPoint(k,ph.GetBinCenter(i),val);g.SetPointError(k,0,float(err));clips+=int(val<0 or val>5);ratios.append({'bin':i,'value':val,'error':float(err)})
  g.SetMarkerStyle(24);g.SetMarkerColor(R.kRed+1);g.SetLineColor(R.kRed+1);g.Draw('P SAME');line=R.TLine(lo,1,hi,1);line.SetLineStyle(2);line.Draw()
  if clips:lab.SetTextSize(.075);lab.DrawLatex(.18,.82,f'{clips} ratio points outside [0,5]')
  c.SaveAs(str(fig/(n+'.pdf')));rf.cd()
  for h in [ph,ih,g]:h.Write()
  rows.append({'variable':n,'covariance_bins':cross.tolist(),'ratios':ratios,'outside_0_5':clips});c.Close()
 rf.Close();pf.Close();fi.Close();subprocess.run(['pdfunite',*[str(fig/(n+'.pdf')) for n in names],str(out/('x_prompt_vs_inclusive_all.pdf' if species=='x' else 'psi2s_prompt_vs_inclusive_all.pdf'))],check=True);report[species]=rows
(out/'prompt_inclusive_comparison.json').write_text(json.dumps({'species':report,'shared_events':len(idx),'errors':'fixed fitted weights and fixed histogram normalizations; overlap covariance sum(w_prompt*w_inclusive)/(norm_prompt*norm_inclusive); no fit/normalization covariance or bootstrap'},indent=2)+'\n')
shutil.copy2(__file__,out/'artifacts'/Path(__file__).name)
print('PROMPT/INCLUSIVE PLOTS PASS',flush=True)
