#!/usr/bin/env python3
from pathlib import Path
import hashlib,json,shutil
import numpy as np
import uproot
import ROOT as R
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetTextFont(42)
repo=Path(__file__).resolve().parents[2];out=repo/'plotER/results/pb18_score_sideband';assert not out.exists();out.mkdir();meta=out/'metadata';meta.mkdir()
old=repo/'plotER/results/pb18_an_score/metadata';data_record=json.loads((old/'Qmumu_manifest.json').read_text());mc_record=json.loads((repo/'plotER/results/prompt_x_pthat5/metadata/provenance.json').read_text())['2018']
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
previous=(old/'plot_pb18_qmumu.py').read_text();cpp=previous.split("R.gInterpreter.Declare('''",1)[1].split("''')",1)[0];R.gInterpreter.Declare(cpp)
inputs={'MC':(Path(mc_record['path']),'ntmix_X3872',mc_record['sha256']),'DATA':(Path(data_record['source']),'ntmix',data_record['source_sha256'])}
edges=np.linspace(-.5,.5,101);histograms={};report={'base_cut':'Q_mumu < 0.13','Q_definition':'double(Bmass) - double(Bmumumass_recovered) - double(Btktkmass)','score_cut':None,'DATA_sideband':'(Bmass > 3.80 && Bmass < 3.85) || (Bmass > 3.89 && Bmass < 3.94)','normalization':'unit area separately','bin_width':.01,'samples':{}}
root=R.TFile.Open(str(meta/'histograms.root'),'RECREATE');keep=[]
for label,(path,tree,expected) in inputs.items():
 before=path.stat();digest=sha(path);assert digest==expected
 columns=['Bmass','Btktkmass','Score']+[f'Bmu{i}{k}' for i in [1,2] for k in ['pt','eta','phi']]+(['pThatreweight'] if label=='MC' else [])
 with uproot.open(path) as f:a=f[tree].arrays(columns,library='np')
 assert all(np.isfinite(z).all() for z in a.values());vec=[]
 for i in [1,2]:
  pt=a[f'Bmu{i}pt'].astype('float64');eta=a[f'Bmu{i}eta'].astype('float64');phi=a[f'Bmu{i}phi'].astype('float64');vec.append(tuple(v.astype('float32') for v in [pt*np.cos(phi),pt*np.sin(phi),pt*np.sinh(eta),np.sqrt((pt*np.cosh(eta))**2+.10565837**2)]))
 sums=[(vec[0][i]+vec[1][i]).astype('float64') for i in range(4)];mass2=sums[3]**2-sum(z*z for z in sums[:3]);assert (mass2>0).all();mm=np.sqrt(mass2).astype('float32');q=a['Bmass'].astype('float64')-mm.astype('float64')-a['Btktkmass'].astype('float64');base=q<.13;mask=base.copy();selection='Q_mumu < 0.13'
 if label=='DATA':
  assert int((base&(a['Score']>.06)).sum())==data_record['selected_entries']
  with uproot.open(old/'Qmumu_recovered.root') as f:assert np.array_equal(q,f['recovered']['Q_mumu'].array(library='np'))
  mask&=((a['Bmass']>3.80)&(a['Bmass']<3.85))|((a['Bmass']>3.89)&(a['Bmass']<3.94));selection+=' && ('+report['DATA_sideband']+')'
 weights=a['pThatreweight'].astype(float) if label=='MC' else np.ones(len(q));assert (weights>=0).all();score=a['Score'][mask].astype(float);w=weights[mask];counts=np.histogram(score,edges,weights=w)[0];squares=np.histogram(score,edges,weights=w*w)[0]
 assert ((score>=edges[0])&(score<edges[-1])).all()
 df=R.RDataFrame(tree,str(path)).Define('Bmumumass_recovered','pb18_restore_mumu(Bmu1pt,Bmu1eta,Bmu1phi,Bmu2pt,Bmu2eta,Bmu2phi)').Define('Q_mumu','double(Bmass)-double(Bmumumass_recovered)-double(Btktkmass)')
 filtered=df.Filter(selection);assert int(filtered.Count().GetValue())==int(mask.sum())
 hptr=filtered.Histo1D((label+'_raw','',100,-.5,.5),'Score','pThatreweight') if label=='MC' else filtered.Histo1D((label+'_raw','',100,-.5,.5),'Score')
 h=hptr.GetValue();assert np.allclose([h.GetBinContent(i) for i in range(1,101)],counts,rtol=1e-8,atol=1e-8);assert np.allclose([h.GetBinError(i)**2 for i in range(1,101)],squares,rtol=1e-8,atol=1e-8);assert h.GetBinContent(0)==h.GetBinContent(101)==0
 root.cd();h.Write();hist=h.Clone(label+'_unit_area');hist.SetDirectory(0);hist.Scale(1/w.sum());hist.Write();histograms[label]=hist;keep.append(hptr)
 report['samples'][label]={'path':str(path),'tree':tree,'sha256':digest,'input_entries':len(q),'base_cut_entries':int(base.sum()),'selected_entries':int(mask.sum()),'score_range':[float(score.min()),float(score.max())],'sumw':float(w.sum()),'sumw2':float(w@w),'weights':'original pThatreweight, including source5' if label=='MC' else 'unit counts','edges':edges.tolist(),'bin_sumw':counts.tolist(),'bin_sumw2':squares.tolist(),'ROOT_numpy_selection_bins_sumw2':'PASS'}
 after=path.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
c=R.TCanvas('score','',900,700);c.SetLeftMargin(.14);c.SetRightMargin(.04);c.SetBottomMargin(.13);c.SetTopMargin(.08)
h=histograms['MC'];h.SetLineColor(R.kBlue+1);h.SetMarkerColor(R.kBlue+1);h.SetMarkerStyle(20);h.SetLineWidth(2);h.SetMarkerSize(.6);h.SetMinimum(0);h.SetMaximum(2.0*max(x.GetMaximum() for x in histograms.values()));h.GetXaxis().SetTitle('Score (BDT)');h.GetYaxis().SetTitle('Unit-area candidates / 0.01');h.GetYaxis().SetTitleOffset(1.5);h.Draw('E1 HIST')
d=histograms['DATA'];d.SetLineColor(R.kOrange+7);d.SetMarkerColor(R.kOrange+7);d.SetLineWidth(2);d.SetMarkerStyle(24);d.SetMarkerSize(.6);d.Draw('E1 HIST SAME')
legend=R.TLegend(.18,.72,.62,.87);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.03);legend.AddEntry(h,'2018 prompt X MC (weighted)','lp');legend.AddEntry(d,'2018 DATA sidebands','lp');legend.Draw()
l=R.TLatex();l.SetNDC();l.SetTextFont(42);l.SetTextSize(.035);l.DrawLatex(.14,.95,'PbPb 2018: prompt X MC vs DATA sidebands');l.SetTextSize(.026)
for y,t in [(.65,'Q = m_{#mu#mu#pi#pi} - m_{#mu#mu} - m_{#pi#pi} < 0.13 GeV/c^{2}'),(.60,'DATA: 3.80 < m < 3.85 or 3.89 < m < 3.94 GeV/c^{2}'),(.55,'No Score cut; each distribution normalized to unit area')]:l.DrawLatex(.18,y,t)
c.RedrawAxis();c.Print(str(out/'Score_MC_DATA_sideband.pdf'));root.Close()
report['status']='PASS';report['error_definition']='sqrt(sumw2)/sumw with fixed normalization';report['script_sha256']=sha(Path(__file__));report['ROOT_version']=R.gROOT.GetVersion();(meta/'manifest.json').write_text(json.dumps(report,indent=2)+'\n');shutil.copy2(__file__,meta/Path(__file__).name)
(out/'README.md').write_text('# 2018 prompt X MC与DATA侧带的Score分布\n\n[Score对比图](Score_MC_DATA_sideband.pdf)。使用Henrique legacy_run2/REflated原样本，MC使用原pThatreweight（保留来源5），DATA单位计数，各自归一化到单位面积。\n\n基础cut沿用最后一版score>0.06质量谱，去掉score cut后为`Bmass - m(mumu)_recovered - Btktkmass < 0.13`；不额外加Bpt、By、中心度、PVz或MC质量窗。DATA另加`3.80 < Bmass < 3.85 or 3.89 < Bmass < 3.94`。\n\n分箱0.01，范围−0.5至0.5，无溢出。误差sqrt(sumw2)/总权重。恢复双muon质量沿用已验收Float_t四动量算法；DATA恢复Q逐条与先前缓存一致，重新加回Score>0.06得到原19456候选。独立C++/ROOT恢复与NumPy的选后计数、所有bin与误差平方一致。[输入、计数及核验](metadata/manifest.json)。\n')
print(json.dumps({k:{n:v for n,v in r.items() if n in ['selected_entries','base_cut_entries','sumw','score_range']} for k,r in report['samples'].items()},indent=2));print('PASS')
