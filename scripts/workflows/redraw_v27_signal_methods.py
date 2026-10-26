#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""User-requested v27 signal/MC and sideband selection comparisons; no refit."""
from pathlib import Path
import json,math,shutil,subprocess,hashlib,tempfile,re
import numpy as np,uproot,numexpr,ROOT as R
from fit_input_contract import selection_expression
from ppref_psi2s_pthat_validation import fingerprint,dump
repo=Path(__file__).resolve().parents[2];base=repo/'plotER/Validation';old=base/'v27_xeff48_signal_methods';oldmc=base/'v27_xeff48_signal_methods_mc'
assert not (old/'archive').exists()
s=json.loads((old/'summary.json').read_text());names=s['variables'];assert len(names)==28 and names[0]=='Prediction'
d=dict(np.load(old/'artifacts/selected_data.npz'));mc=dict(np.load(oldmc/'artifacts/mc_selected.npz'))
prov=json.loads((old/'provenance.json').read_text());source=Path(prov['data']['path']);assert fingerprint(source)['sha256']==prov['data']['sha256']
cols=sorted(set(names+['Bmass']+re.findall(r'\bB\w+\b',s['selection'])))
with uproot.open(source) as f:a=f['ntmix'].arrays(cols,library='np')
mask=numexpr.evaluate(selection_expression(s['selection']),local_dict=a)&(((a['Bmass']>3.75)&(a['Bmass']<3.84))|((a['Bmass']>3.90)&(a['Bmass']<3.99)))
pre={k:v[mask] for k,v in a.items()};pre['source_entry']=np.flatnonzero(mask)
post=pre['Prediction']>s['threshold'];assert np.array_equal(pre['source_entry'][post],d['source_entry'][d['sideband_weight']>0])
for n in names:assert np.array_equal(pre[n][post],d[n][d['sideband_weight']>0])
protected={label:{str(p.relative_to(root)):fingerprint(p)['sha256'] for p in root.rglob('*') if p.is_file()} for label,root in [('three_curve',old),('four_curve',oldmc)]}
work=Path(tempfile.mkdtemp(prefix='v27_redraw_',dir='/tmp/leyao'));print('WORK',work,flush=True)
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0)
R.gInterpreter.Declare('#include "'+str(repo/'plotER/validation_archive/aux.h')+'"');titles={str(q.expr):str(q.title).split(';')[1] for q in R.getSignalVars('ntmix_X3872')};titles['Prediction']='Prediction (XGBoost score)'
original_metrics=json.loads((old/'variable_metrics.json').read_text());ranges={r['variable']:(*r['range'],r['absolute_display']) for r in original_metrics}
for study in ['signal_mc','sideband_selection']:
 out=work/study;out.mkdir();art=out/'artifacts';art.mkdir();(out/'pdf').mkdir()
 if study=='signal_mc':
  weights={'splot':d['splot_weight'],'subtraction':d['subtraction_weight'],'mc':mc['pThatreweight'].astype(float)};samples={'splot':d,'subtraction':d,'mc':mc};denominator='mc';ratio_title='Method / MC'
 else:
  weights={'splot':d['splot_weight'],'sideband':d['sideband_weight'],'pre_sideband':np.ones(len(pre['Prediction']))};samples={'splot':d,'sideband':d,'pre_sideband':pre};denominator='pre_sideband';ratio_title='Method / pre-ML SB';ranges['Prediction']=(0.,1.,False)
 root=R.TFile(str(art/'histograms.root'),'RECREATE');records=[];colors={'splot':R.kRed+1,'subtraction':R.kGreen+2,'sideband':R.kBlue+1,'mc':R.kOrange+7,'pre_sideband':R.kViolet+1};labels={'splot':'X sPlot','subtraction':'Signal - sideband/3','sideband':'Sideband, score cut','mc':'Prompt X MC (pThat)','pre_sideband':'Sideband, pre-ML cuts'}
 for j,n in enumerate(names):
  lo,hi,ab=ranges[n];x=abs(d[n]) if ab else d[n];edges=np.linspace(lo,hi,16);inside=(x>=lo)&(x<hi);hs={};norms={};raws={};rec=dict(variable=n,range=[lo,hi],absolute_display=ab,samples={})
  for k,weight in weights.items():
   vals=samples[k][n];x=abs(vals) if ab else vals;inside=(x>=lo)&(x<hi)
   binidx=np.searchsorted(edges,x[inside],side='right')-1;raw=np.bincount(binidx,weights=weight[inside],minlength=15);err2=np.bincount(binidx,weights=weight[inside]**2,minlength=15);norm=raw.sum();fallback=norm<=0
   if fallback:norm=weight.sum()
   assert norm>0;norms[k]=norm;raws[k]=raw;h=R.TH1D(n+'_'+k,'',15,lo,hi);h.SetDirectory(0);h.Sumw2()
   for i in range(15):h.SetBinContent(i+1,raw[i]/norm);h.SetBinError(i+1,math.sqrt(err2[i])/norm)
   h.SetLineColor(colors[k]);h.SetMarkerColor(colors[k]);h.SetMarkerStyle({'splot':24,'subtraction':25,'sideband':20,'mc':1,'pre_sideband':26}[k]);h.SetLineWidth(2);hs[k]=h;rec['samples'][k]=dict(normalization=float(norm),fallback=bool(fallback),inrange_fraction=float(raw.sum()/weight.sum()),raw_bins=raw.tolist(),raw_variances=err2.tolist())
  c=R.TCanvas('c'+str(j),'',760,650);top=R.TPad('top'+str(j),'',0,.30,1,1);bot=R.TPad('bot'+str(j),'',0,0,1,.30)
  for pad in [top,bot]:pad.SetLeftMargin(.14);pad.SetRightMargin(.04);pad.Draw()
  top.SetBottomMargin(.01);bot.SetTopMargin(.01);bot.SetBottomMargin(.33);top.cd();h=hs['splot'];high=max(hh.GetBinContent(i)+hh.GetBinError(i) for hh in hs.values() for i in range(1,16));low=min(0,min(hh.GetBinContent(i)-hh.GetBinError(i) for hh in hs.values() for i in range(1,16)));h.SetMaximum(max(high,.1)*1.8);h.SetMinimum(low*1.15);h.GetYaxis().SetTitle('Normalized entries');h.GetXaxis().SetLabelSize(0);h.Draw('E')
  for k,hh in hs.items():
   if k=='mc':
    band=hh.Clone(n+'_mc_band');band.SetDirectory(0);band.SetFillColorAlpha(R.kOrange+7,.25);band.SetMarkerSize(0);band.Draw('E2 SAME');hh.Draw('HIST SAME')
   else:hh.Draw('E SAME')
  leg=R.TLegend(.57,.70,.95,.92);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.033)
  for k in weights:leg.AddEntry(hs[k],labels[k],'l' if k=='mc' else 'lep')
  leg.Draw();lab=R.TLatex();lab.SetNDC();lab.SetTextSize(.047);lab.DrawLatex(.18,.87,'PbPb23 X(3872)');lab.SetTextSize(.030);lab.DrawLatex(.18,.80,'v27 v1, optimal X efficiency 48%');lab.DrawLatex(.18,.74,'Signal: (3.84,3.90)');lab.DrawLatex(.18,.68,'SB: (3.75,3.84) + (3.90,3.99)')
  bot.cd();axis=R.TH1D(n+'_axis','',15,lo,hi);axis.SetDirectory(0);axis.SetMinimum(0);axis.SetMaximum(5);axis.GetYaxis().SetTitle(ratio_title);axis.GetYaxis().SetTitleSize(.085);axis.GetYaxis().SetLabelSize(.08);axis.GetYaxis().SetTitleOffset(.75);axis.GetYaxis().SetNdivisions(304);axis.GetXaxis().SetTitle(titles.get(n,n));axis.GetXaxis().SetTitleSize(.11);axis.GetXaxis().SetLabelSize(.10);axis.Draw('AXIS');graphs=[];clips=0
  x=abs(d[n]) if ab else d[n];inside=(x>=lo)&(x<hi)
  for k in [key for key in weights if key!=denominator]:
   g=R.TGraphErrors();g.SetName(n+'_'+k+'_ratio');g.SetMarkerStyle(24 if k=='splot' else 25);g.SetMarkerColor(colors[k]);g.SetLineColor(colors[k]);cov=np.zeros(15) if denominator=='mc' else np.histogram(x[inside],edges,weights=(weights[k]*d['sideband_weight'])[inside])[0]/(norms[k]*norms[denominator])
   for i in range(15):
    den=hs[denominator].GetBinContent(i+1);num=hs[k].GetBinContent(i+1)
    if den<=0:continue
    val=num/den;var=hs[k].GetBinError(i+1)**2/den**2+num*num*hs[denominator].GetBinError(i+1)**2/den**4-2*num*cov[i]/den**3
    if denominator=='pre_sideband':
     inbin=(x>=edges[i])&(x<edges[i+1]);overlap=d['sideband_weight'][inbin];residual=weights[k][inbin]/norms[k]-val*overlap/norms[denominator];extra=raws[denominator][i]-overlap.sum();assert extra>=0;var=(np.square(residual).sum()+extra*(val/norms[denominator])**2)/den**2
    assert var>=0;idx=g.GetN();g.SetPoint(idx,hs[k].GetBinCenter(i+1),val);g.SetPointError(idx,0,math.sqrt(max(0,var)));clips+=int(val<0 or val>5)
   g.Draw('P SAME');graphs.append(g)
  line=R.TLine(lo,1,hi,1);line.SetLineStyle(2);line.Draw()
  if clips:lab.SetTextSize(.07);lab.DrawLatex(.18,.85,f'{clips} ratio points outside [0,5]')
  rec['ratio_points_outside_display']=clips;c.SaveAs(str(out/'pdf'/f'{n}.pdf'));root.cd()
  for obj in [*hs.values(),*graphs]:obj.Write()
  c.Close();records.append(rec)
 root.Close();dump(out/'variable_metrics.json',records);subprocess.run(['pdfunite',*[str(out/'pdf'/f'{n}.pdf') for n in names],str(work/(study+'_all.pdf'))],check=True)
 dump(out/'summary.json',{**s,'study':study,'samples':list(weights),'denominator':denominator,'prediction_range':ranges['Prediction'][:2],'pre_sideband_entries':len(pre['Prediction']),'post_sideband_entries':int(post.sum()),'refit':False,'source_data':prov['data'],'mc_weight':'pThatreweight','mc_pthat5_retained':True,'mc_mass_cut':None})
 np.savez_compressed(out/'artifacts/pre_sideband.npz',**pre) if study=='sideband_selection' else None
 shutil.copy2(Path(__file__),art/Path(__file__).name)
for label,root in [('three_curve',old),('four_curve',oldmc)]:
 for rel,h in protected[label].items():assert fingerprint(root/rel)['sha256']==h
archive=old/'archive';archive.mkdir();dest=archive/'splot_sbs_sideband';dest.mkdir()
for p in list(old.iterdir()):
 if p.name!='archive':shutil.move(str(p),dest/p.name)
shutil.move(str(oldmc),archive/'splot_sbs_mc_sideband')
for p in list(work.iterdir()):shutil.move(str(p),old/p.name)
work.rmdir();dump(archive/'original_hashes.json',protected)
for label,folder in [('three_curve','splot_sbs_sideband'),('four_curve','splot_sbs_mc_sideband')]:
 for rel,h in protected[label].items():assert fingerprint(archive/folder/rel)['sha256']==h
print('COMPLETE',len(pre['Prediction']),int(post.sum()),flush=True)
