#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""v27 xeff48 X sPlot versus width-scaled sideband subtraction and raw sideband."""
from pathlib import Path
import sys,json,hashlib,math,re,shutil,subprocess
import numpy as np,uproot,numexpr,ROOT as R
from fit_input_contract import selection_expression
from ppref_psi2s_pthat_validation import fingerprint,dump,record,stats
repo=Path(__file__).resolve().parents[2];tag='X_pb23_v27_fid18_9v9_rw0_xgb_v1';source=repo/'fitER/results/ml_fits'/tag/'two_peak_mc_width/xeff48';out=repo/'plotER/Validation/v27_xeff48_signal_methods'
assert not out.exists()
names='Bmu1pt Bmu2pt Btrk1Pt Btrk2Pt Btrk1dR Btrk2dR BtrkPtimb Bmu1eta Bmu2eta Btrk1Eta Btrk2Eta Bujeta Btktkpt Btktketa Bchi2Prob BujvProb BtktkvProb Balpha Bcos_dtheta Bnorm_svpvDistance_2D BsvpvDistance_2D Bdtheta Bujpt Btrk1Dxy1 Btrk2Dxy1 Bnorm_trk1Dxy Bnorm_trk2Dxy'.split();assert len(names)==27
m=json.loads((source/'manifest.json').read_text());ref=json.loads((source/'fit_result.json').read_text());assert ref['audit']['quality_pass']
paths={'data':Path(m['inputs']['data']['source']),'workspace':source/'fit_workspace.root','manifest':source/'manifest.json','fit_result':source/'fit_result.json','script':Path(__file__)};provenance={k:fingerprint(p) for k,p in paths.items()}
cols=sorted(set(names+['Bmass','Prediction']+re.findall(r'\bB\w+\b',m['selection'])))
with uproot.open(paths['data']) as f:a=f['ntmix'].arrays(cols,library='np')
mask=numexpr.evaluate(selection_expression(m['selection']),local_dict=a)&(a['Prediction']>m['threshold'])&(a['Bmass']>3.62)&(a['Bmass']<4)
d={n:a[n][mask] for n in cols};d['source_entry']=np.flatnonzero(mask);assert len(d['Bmass'])==ref['data_entries'];assert all(np.isfinite(v).all() for v in d.values())
cache=source/'artifacts/data_cache.root';assert fingerprint(cache)['sha256']==m['inputs']['data']['sha256']
with uproot.open(cache) as f:
 for n,v in f['ntmix'].arrays(library='np').items():assert np.array_equal(d[n],v)
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0);assert R.gROOT.GetVersion()=='6.32.02';R.RooMsgService.instance().setGlobalKillBelow(R.RooFit.WARNING)
out.mkdir();art=out/'artifacts';art.mkdir();(out/'pdf').mkdir();dump(out/'provenance.json',provenance)
f=R.TFile.Open(str(paths['workspace']));w=f.Get('ws_nominal');model=w.pdf('model');data=w.data('data_data');mass=w.var('Bmass');pars=model.getParameters(data);pars.assignValueOnly(f.Get('fit_result_alt').floatParsFinal())
for i,x in enumerate(d['Bmass']):assert data.get(i).getRealValue('Bmass')==x
ykeys=['x_yield','psi_yield','nbkg'];yields=[w.var(k) for k in ykeys]
for p in pars:p.setConstant(p.GetName() not in ykeys)
fixed={p.GetName():p.getVal() for p in pars if p.isConstant()};assert fixed['x_scale']==fixed['psi_scale']==1
opts=[R.RooFit.Save(),R.RooFit.Extended(True),R.RooFit.Range('all'),R.RooFit.Strategy(2),R.RooFit.Hesse(True),R.RooFit.PrintLevel(-1)]
good=lambda q:q.status()==0 and q.covQual()==3 and math.isfinite(q.edm()) and q.edm()<1e-3
fit=model.fitTo(data,*opts);history=[record(fit)]
if not good(fit) and math.isfinite(fit.minNll()):pars.assignValueOnly(fit.floatParsFinal());fit=model.fitTo(data,*opts);history.append(record(fit))
dump(out/'retry_history.json',history);assert good(fit);assert set(p.GetName() for p in fit.floatParsFinal())==set(ykeys);assert all(not p['near_boundary'] for p in record(fit)['parameters'].values())
pre={p.GetName():p.getVal() for p in yields};original_min={p.GetName():p.getMin() for p in yields}
for p in yields:p.setMin(0.)
spl=R.RooStats.SPlot('splot','splot',data,model,R.RooArgList(*yields),R.RooArgSet(),True,False,'',R.RooFit.Strategy(2),R.RooFit.Hesse(True),R.RooFit.PrintLevel(-1),R.RooFit.Range('all'))
sw=np.array([[spl.GetSWeight(i,k) for k in ykeys] for i in range(data.numEntries())]);assert np.isfinite(sw).all()
for p in yields:
 assert p.getVal()>original_min[p.GetName()];p.setMin(original_min[p.GetName()]);assert abs(p.getVal()-pre[p.GetName()])<.01*fit.floatParsFinal().find(p.GetName()).getError()
for n,v in fixed.items():assert w.var(n).isConstant() and w.var(n).getVal()==v
components=[w.pdf(k) for k in ['x_signal','psi_signal','backgroundPdf']];obs=R.RooArgSet(mass);P=[]
for x in d['Bmass']:mass.setVal(float(x));P.append([p.getVal(obs) for p in components])
P=np.array(P);nv=np.array([p.getVal() for p in yields]);den=P@nv;A=P/den[:,None];V=np.linalg.inv(A.T@A);assert np.max(abs(sw-A@V.T))<1e-8;assert np.max(abs(sw.sum(axis=0)-nv))<1e-5*max(nv);assert np.max(abs(sw.sum(axis=1)-1))<1e-3
massvalues=d['Bmass'];sig=(massvalues>3.84)&(massvalues<3.90);left=(massvalues>3.75)&(massvalues<3.84);right=(massvalues>3.90)&(massvalues<3.99);sb=left|right;assert not np.any(sig&sb)
weights={'splot':sw[:,0],'subtraction':sig.astype(float)-sb.astype(float)/3,'sideband':sb.astype(float)}
for k,v in weights.items():assert v.sum()>0;d[k+'_weight']=v
np.savez_compressed(art/'selected_data.npz',**d);np.savez(art/'sweight_algebra.npz',component_pdf=P,yields=nv,covariance=V,sweights=sw)
rf=R.TFile(str(art/'splot_workspace.root'),'RECREATE');w.Write();data.Write('splot_data');fit.Write('fit_yields');rf.Close()
summary=dict(tag=tag,xeff=48,threshold=m['threshold'],selection=m['selection'],variables=names,fit_source=str(source),fit=record(fit),fixed_parameters=fixed,splot_yields=dict(zip(ykeys,nv.tolist())),counts=dict(data=len(sig),signal=int(sig.sum()),left=int(left.sum()),right=int(right.sum())),background_width_ratio=1/3,signal_window=[3.84,3.90],sideband_windows=[[3.75,3.84],[3.90,3.99]],weight_stats={k:stats(v) for k,v in weights.items()},checks=dict(source_entry_alignment=True,analytic_sweight_max_difference=float(np.max(abs(sw-A@V.T))),event_sum_max_difference=float(np.max(abs(sw.sum(axis=1)-1)))),limitations=['sPlot assumes mass-variable independence within components after score cut.','Subtraction assumes sideband variable shape represents signal-window background; no signal-tail correction.','Negative bins and signed weights retained.','Histogram errors use sumw2 with fixed fitted weights and fixed normalization; no shape/normalization uncertainty or bootstrap.','Ratio errors include same-event numerator/sideband covariance with fixed weights.'])
dump(out/'summary.json',summary)
template=json.loads((repo/'plotER/validation_archive/archive/previous_validation/results/pb23_x_mc_sideband/variable_metrics.json').read_text());ranges={r['variable']:(*r['range'],r['absolute_display']) for r in template}
for rel in ['plotER/Validation/pb23_x_mc_sideband_no_dr/artifacts/vertex_tracks/variable_metrics.json']:
 for r in json.loads((repo/rel).read_text()):ranges[r['variable']]=(*r['range'],r['absolute_display'])
R.gInterpreter.Declare('#include "'+str(repo/'plotER/validation_archive/aux.h')+'"');titles={str(q.expr):str(q.title).split(';')[1] for q in R.getSignalVars('ntmix_X3872')}
root=R.TFile(str(art/'histograms.root'),'RECREATE');records=[];colors={'splot':R.kRed+1,'subtraction':R.kGreen+2,'sideband':R.kBlue+1};labels={'splot':'X sPlot','subtraction':'Signal - sideband/3','sideband':'Sideband DATA'}
for j,n in enumerate(names):
 lo,hi,ab=ranges[n];x=abs(d[n]) if ab else d[n];edges=np.linspace(lo,hi,16);inside=(x>=lo)&(x<hi);hs={};norms={};raws={};rec=dict(variable=n,range=[lo,hi],absolute_display=ab,samples={})
 for k,weight in weights.items():
  raw=np.histogram(x[inside],edges,weights=weight[inside])[0];err2=np.histogram(x[inside],edges,weights=weight[inside]**2)[0];norm=raw.sum();fallback=norm<=0
  if fallback:norm=weight.sum()
  assert norm>0;norms[k]=norm;raws[k]=raw;h=R.TH1D(n+'_'+k,'',15,lo,hi);h.SetDirectory(0);h.Sumw2()
  for i in range(15):h.SetBinContent(i+1,raw[i]/norm);h.SetBinError(i+1,math.sqrt(err2[i])/norm)
  h.SetLineColor(colors[k]);h.SetMarkerColor(colors[k]);h.SetMarkerStyle({'splot':24,'subtraction':25,'sideband':20}[k]);h.SetLineWidth(2);hs[k]=h;rec['samples'][k]=dict(normalization=float(norm),fallback=bool(fallback),inrange_fraction=float(raw.sum()/weight.sum()),raw_bins=raw.tolist(),raw_variances=err2.tolist())
 c=R.TCanvas('c'+str(j),'',760,650);top=R.TPad('top'+str(j),'',0,.30,1,1);bot=R.TPad('bot'+str(j),'',0,0,1,.30)
 for pad in [top,bot]:pad.SetLeftMargin(.14);pad.SetRightMargin(.04);pad.Draw()
 top.SetBottomMargin(.01);bot.SetTopMargin(.01);bot.SetBottomMargin(.33);top.cd();h=hs['splot'];high=max(hh.GetBinContent(i)+hh.GetBinError(i) for hh in hs.values() for i in range(1,16));low=min(0,min(hh.GetBinContent(i)-hh.GetBinError(i) for hh in hs.values() for i in range(1,16)));h.SetMaximum(max(high,.1)*1.8);h.SetMinimum(low*1.15);h.GetYaxis().SetTitle('Normalized entries');h.GetXaxis().SetLabelSize(0);h.Draw('E');hs['subtraction'].Draw('E SAME');hs['sideband'].Draw('E SAME')
 leg=R.TLegend(.57,.74,.95,.92);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.033)
 for k in weights:leg.AddEntry(hs[k],labels[k],'lep')
 leg.Draw();lab=R.TLatex();lab.SetNDC();lab.SetTextSize(.047);lab.DrawLatex(.18,.87,'PbPb23 X(3872)');lab.SetTextSize(.030);lab.DrawLatex(.18,.80,'v27 v1, X efficiency 48%');lab.DrawLatex(.18,.74,'Signal: (3.84,3.90)');lab.DrawLatex(.18,.68,'SB: (3.75,3.84) + (3.90,3.99)')
 bot.cd();axis=R.TH1D(n+'_axis','',15,lo,hi);axis.SetDirectory(0);axis.SetMinimum(0);axis.SetMaximum(5);axis.GetYaxis().SetTitle('Method / sideband');axis.GetYaxis().SetTitleSize(.085);axis.GetYaxis().SetLabelSize(.08);axis.GetYaxis().SetTitleOffset(.75);axis.GetYaxis().SetNdivisions(304);axis.GetXaxis().SetTitle(titles.get(n,n));axis.GetXaxis().SetTitleSize(.11);axis.GetXaxis().SetLabelSize(.10);axis.Draw('AXIS');graphs=[];clips=0
 for k in ['splot','subtraction']:
  g=R.TGraphErrors();g.SetName(n+'_'+k+'_ratio');g.SetMarkerStyle(24 if k=='splot' else 25);g.SetMarkerColor(colors[k]);g.SetLineColor(colors[k]);cov=np.histogram(x[inside],edges,weights=(weights[k]*weights['sideband'])[inside])[0]/(norms[k]*norms['sideband'])
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
for k,p in paths.items():assert fingerprint(p)==provenance[k]
for p in [Path(__file__),repo/'scripts/workflows/ppref_psi2s_pthat_validation.py',repo/'plotER/validation_archive/aux.h']:shutil.copy2(p,art/p.name)
(out/'README.md').write_text('# v27 xeff48：X sPlot、sideband subtraction与sideband\n\n[27页合订本](validation_all.pdf) · [单图](pdf/) · [计数与sPlot验收](summary.json)\n\n共同选择：'+m['selection']+f"，Prediction>{m['threshold']:.12g}。sPlot来自3.62–4.00双峰，冻结全部峰形和背景形状后浮动三个yield。\n\nSignal为(3.84,3.90)，sideband为(3.75,3.84)/(3.90,3.99)，严格开区间；subtraction逐bin为Nsignal−Nsideband/3，统计方差为Nsignal+Nsideband/9。三条曲线各自按可见范围积分归一化，非正积分时使用全范围净权重；保留负权重/负bin。\n\n误差按固定权重SumW2；比值包括同一事例引起的协方差，但不包含拟合形状与归一化常数的不确定性。sPlot依赖组分内质量与变量独立，subtraction依赖侧带能代表信号窗背景，未校正信号尾进入侧带。两种提取使用相同DATA，并非独立测量。无MC曲线、bootstrap或额外拟合扫描。\n")
print('COMPLETE',json.dumps(summary['counts']),json.dumps(summary['weight_stats']),flush=True)
