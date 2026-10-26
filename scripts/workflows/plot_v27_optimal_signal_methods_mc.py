#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""Add prompt MC to accepted v27 xeff48 signal-method plots without refitting."""
from pathlib import Path
import json,math,re,shutil,subprocess,hashlib
import numpy as np,uproot,numexpr,ROOT as R
from fit_input_contract import load_training,selection_expression
from ppref_psi2s_pthat_validation import fingerprint,dump,stats
repo=Path(__file__).resolve().parents[2];old=repo/'plotER/Validation/v27_xeff48_signal_methods';out=repo/'plotER/Validation/v27_xeff48_signal_methods_mc';assert not out.exists()
protected={str(p.relative_to(old)):fingerprint(p)['sha256'] for p in old.rglob('*') if p.is_file()}
s=json.loads((old/'summary.json').read_text());names=s['variables'];d=dict(np.load(old/'artifacts/selected_data.npz'));weights={k:d[k+'_weight'] for k in ['splot','subtraction','sideband']}
t,points=load_training(repo,s['tag']);assert t['selection']==s['selection'];src=Path(t['sources']['mc']['path']);prov=fingerprint(src)
cols=sorted(set(names+['Bmass','Prediction','pThatreweight']+re.findall(r'\bB\w+\b',t['selection'])))
with uproot.open(src) as f:a=f[t['sources']['mc']['tree']].arrays(cols,library='np')
mask=numexpr.evaluate(selection_expression(t.get('mc_apply_selection',t['selection'])),local_dict=a)&(a['Prediction']>s['threshold'])
mc={k:v[mask] for k,v in a.items()};mc['source_entry']=np.flatnonzero(mask);assert all(np.isfinite(v).all() for v in mc.values());weights['mc']=mc['pThatreweight'].astype(float);assert (weights['mc']>0).all();assert np.any(abs(weights['mc']-1.8475149726206412)<1e-12)
out.mkdir();art=out/'artifacts';art.mkdir();(out/'pdf').mkdir();np.savez_compressed(art/'mc_selected.npz',**mc);dump(art/'original_hashes.json',protected)
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0)
template=json.loads((repo/'plotER/validation_archive/archive/previous_validation/results/pb23_x_mc_sideband/variable_metrics.json').read_text());ranges={r['variable']:(*r['range'],r['absolute_display']) for r in template}
for rel in ['plotER/Validation/pb23_x_mc_sideband_no_dr/artifacts/vertex_tracks/variable_metrics.json']:
 for r in json.loads((repo/rel).read_text()):ranges[r['variable']]=(*r['range'],r['absolute_display'])
R.gInterpreter.Declare('#include "'+str(repo/'plotER/validation_archive/aux.h')+'"');titles={str(q.expr):str(q.title).split(';')[1] for q in R.getSignalVars('ntmix_X3872')}
root=R.TFile(str(art/'histograms.root'),'RECREATE');records=[];colors={'splot':R.kRed+1,'subtraction':R.kGreen+2,'sideband':R.kBlue+1,'mc':R.kOrange+7};labels={'splot':'X sPlot','subtraction':'Signal - sideband/3','sideband':'Sideband DATA','mc':'Prompt X MC (pThat)'}
for j,n in enumerate(names):
 lo,hi,ab=ranges[n];x=abs(d[n]) if ab else d[n];edges=np.linspace(lo,hi,16);inside=(x>=lo)&(x<hi);hs={};norms={};raws={};rec=dict(variable=n,range=[lo,hi],absolute_display=ab,samples={})
 for k,weight in weights.items():
  vals=mc[n] if k=='mc' else d[n];x=abs(vals) if ab else vals;inside=(x>=lo)&(x<hi)
  raw=np.histogram(x[inside],edges,weights=weight[inside])[0];err2=np.histogram(x[inside],edges,weights=weight[inside]**2)[0];norm=raw.sum();fallback=norm<=0
  if fallback:norm=weight.sum()
  assert norm>0;norms[k]=norm;raws[k]=raw;h=R.TH1D(n+'_'+k,'',15,lo,hi);h.SetDirectory(0);h.Sumw2()
  for i in range(15):h.SetBinContent(i+1,raw[i]/norm);h.SetBinError(i+1,math.sqrt(err2[i])/norm)
  h.SetLineColor(colors[k]);h.SetMarkerColor(colors[k]);h.SetMarkerStyle({'splot':24,'subtraction':25,'sideband':20,'mc':1}[k]);h.SetLineWidth(2);hs[k]=h;rec['samples'][k]=dict(normalization=float(norm),fallback=bool(fallback),inrange_fraction=float(raw.sum()/weight.sum()),raw_bins=raw.tolist(),raw_variances=err2.tolist())
 c=R.TCanvas('c'+str(j),'',760,650);top=R.TPad('top'+str(j),'',0,.30,1,1);bot=R.TPad('bot'+str(j),'',0,0,1,.30)
 for pad in [top,bot]:pad.SetLeftMargin(.14);pad.SetRightMargin(.04);pad.Draw()
 top.SetBottomMargin(.01);bot.SetTopMargin(.01);bot.SetBottomMargin(.33);top.cd();h=hs['splot'];high=max(hh.GetBinContent(i)+hh.GetBinError(i) for hh in hs.values() for i in range(1,16));low=min(0,min(hh.GetBinContent(i)-hh.GetBinError(i) for hh in hs.values() for i in range(1,16)));h.SetMaximum(max(high,.1)*1.8);h.SetMinimum(low*1.15);h.GetYaxis().SetTitle('Normalized entries');h.GetXaxis().SetLabelSize(0);h.Draw('E');hs['subtraction'].Draw('E SAME');hs['sideband'].Draw('E SAME');band=hs['mc'].Clone(n+'_mc_band');band.SetDirectory(0);band.SetFillColorAlpha(R.kOrange+7,.25);band.SetMarkerSize(0);band.Draw('E2 SAME');hs['mc'].Draw('HIST SAME');hs['splot'].Draw('E SAME');hs['subtraction'].Draw('E SAME');hs['sideband'].Draw('E SAME')
 leg=R.TLegend(.57,.70,.95,.92);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.033)
 for k in ['subtraction','splot','mc','sideband']:leg.AddEntry(hs[k],labels[k],'l' if k=='mc' else 'lep')
 leg.Draw();lab=R.TLatex();lab.SetNDC();lab.SetTextSize(.047);lab.DrawLatex(.18,.87,'PbPb23 X(3872)');lab.SetTextSize(.030);lab.DrawLatex(.18,.80,'v27 v1, X efficiency 48%');lab.DrawLatex(.18,.74,'Signal: (3.84,3.90)');lab.DrawLatex(.18,.68,'SB: (3.75,3.84) + (3.90,3.99)')
 bot.cd();axis=R.TH1D(n+'_axis','',15,lo,hi);axis.SetDirectory(0);axis.SetMinimum(0);axis.SetMaximum(5);axis.GetYaxis().SetTitle('Method / sideband');axis.GetYaxis().SetTitleSize(.085);axis.GetYaxis().SetLabelSize(.08);axis.GetYaxis().SetTitleOffset(.75);axis.GetYaxis().SetNdivisions(304);axis.GetXaxis().SetTitle(titles.get(n,n));axis.GetXaxis().SetTitleSize(.11);axis.GetXaxis().SetLabelSize(.10);axis.Draw('AXIS');graphs=[];clips=0
 x=abs(d[n]) if ab else d[n];inside=(x>=lo)&(x<hi)
 for k in ['splot','subtraction','mc']:
  g=R.TGraphErrors();g.SetName(n+'_'+k+'_ratio');g.SetMarkerStyle(24 if k=='splot' else 25);g.SetMarkerColor(colors[k]);g.SetLineColor(colors[k]);cov=np.zeros(15) if k=='mc' else np.histogram(x[inside],edges,weights=(weights[k]*weights['sideband'])[inside])[0]/(norms[k]*norms['sideband'])
  for i in range(15):
   den=hs['sideband'].GetBinContent(i+1);num=hs[k].GetBinContent(i+1)
   if den<=0:continue
   val=num/den;var=hs[k].GetBinError(i+1)**2/den**2+num*num*hs['sideband'].GetBinError(i+1)**2/den**4-2*num*cov[i]/den**3;assert var>-1e-9;idx=g.GetN();g.SetPoint(idx,hs[k].GetBinCenter(i+1),val);g.SetPointError(idx,0,math.sqrt(max(0,var)));clips+=int(val<0 or val>5)
  g.Draw('P SAME');graphs.append(g)
 line=R.TLine(lo,1,hi,1);line.SetLineStyle(2);line.Draw()
 if clips:lab.SetTextSize(.07);lab.DrawLatex(.18,.85,f'{clips} ratio points outside [0,5]')
 rec['ratio_points_outside_display']=clips;c.SaveAs(str(out/'pdf'/f'{n}.pdf'));root.cd()
 for obj in [*hs.values(),*graphs]:obj.Write()
 c.Close();records.append(rec)
root.Close();dump(out/'variable_metrics.json',records);subprocess.run(['pdfunite',*[str(out/'pdf'/f'{n}.pdf') for n in names],str(out/'validation_all.pdf')],check=True)

assert fingerprint(src)==prov
for rel,h in protected.items():assert fingerprint(old/rel)['sha256']==h
summary={**s,'mc_source':prov,'mc_selection':t.get('mc_apply_selection',t['selection']),'mc_mass_cut':None,'mc_weight':'pThatreweight','mc_pthat5_retained':True,'mc':stats(weights['mc']),'source_three_curve_study':str(old),'old_files_unchanged':True,'refit':False}
dump(out/'summary.json',summary);shutil.copy2(Path(__file__),art/Path(__file__).name)
(out/'README.md').write_text('# v27 optimal48%：SBS、sPlot、prompt MC、SB\n\n[27页合订本](validation_all.pdf) · [单图](pdf/) · [筛选与计数](summary.json)\n\n复用原sPlot/SBS/SB权重与27变量范围、分箱和布局，未重拟合。原三曲线study所有文件保持不变。MC使用相同fid18和Prediction>0.793513834476471，pThatreweight加权、保留pthat5、无额外质量或位移筛选；prompt MC与inclusive DATA信号比较。\n\nSignal (3.84,3.90)，SB (3.75,3.84)/(3.90,3.99)，SBS=signal−SB/3。各曲线按可见范围各自归一化，保留负bin；MC橙色线与统计误差带。下方继续以SB为分母；DATA方法比值含共享事例协方差，MC与DATA视为独立。固定权重/归一化误差不含拟合形状与归一化波动。\n')
print('COMPLETE MC',json.dumps(summary['mc']))
