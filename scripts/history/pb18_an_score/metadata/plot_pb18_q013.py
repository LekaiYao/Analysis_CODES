from pathlib import Path
import json,hashlib,shutil,array
import numpy as np,uproot,ROOT as R
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0)
repo=Path('/eos/home-l/leyao/pbpb_work/X_analysis/Analysis_CODES');out=Path('/tmp/pb18_q013_review');assert not out.exists();(out/'metadata').mkdir(parents=True)
src=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_DATA.root');st=src.stat();ref=repo/'docs/references/AN2019_021_v10.pdf'
with uproot.open(src) as f:a=f['ntmix'].arrays(['Bmass','Score','BQvalue'],library='np')
assert all(np.isfinite(x).all() for x in a.values());mask=(a['Score']>.06)&(a['BQvalue']<.13);edges=np.linspace(3.6,4.0,81);count,_=np.histogram(a['Bmass'][mask],bins=edges)
# ROOT defines its upper histogram edge as overflow.
count[-1]-=int(np.sum(mask&(a['Bmass']==4.0)))
f=R.TFile.Open(str(src));t=f.Get('ntmix');R.gROOT.cd();h=R.TH1D('Bmass_score0p06',';m_{J/#psi #pi^{+}#pi^{-}} [GeV/c^{2}];Entries / 5 MeV/c^{2}',80,3.6,4.0);h.Sumw2();n=t.Draw('Bmass>>Bmass_score0p06','Score > 0.06 && BQvalue < 0.13','goff');assert n==int(mask.sum());assert np.array_equal(count,np.array([int(h.GetBinContent(i)) for i in range(1,81)]));assert h.GetBinContent(0)==int(np.sum(mask&(a['Bmass']<3.6)));assert h.GetBinContent(81)==int(np.sum(mask&(a['Bmass']>=4.0)))
c=R.TCanvas('cPb18AnScore','',1000,650);c.SetLeftMargin(.10);c.SetRightMargin(.36);c.SetTopMargin(.06);c.SetBottomMargin(.12)
h.SetLineColor(R.kBlue+1);h.SetLineWidth(2);h.SetFillColor(R.kBlue);h.SetFillStyle(3358);h.SetMinimum(0);h.SetMaximum(1.25*h.GetMaximum());h.Draw('HIST')
l=R.TLegend(.66,.55,.985,.91);l.SetBorderSize(0);l.SetFillColor(R.kWhite);l.SetTextSize(.032);l.SetHeader('#bf{PbPb 2018 Data}');l.AddEntry(h,f'N = {n:,}','f')
for label in ['Existing REflated sample','Score (BDT) > 0.06','AN-19-021, Table 3','BQvalue < 0.13 GeV/c^{2}'] :l.AddEntry(R.nullptr,label,'')
l.Draw();c.SaveAs(str(out/'Bmass_score0p06.pdf'))
r=R.TFile.Open(str(out/'metadata/Bmass_score0p06.root'),'RECREATE');h.Write();r.Close();f.Close()
np.savetxt(out/'metadata/Bmass_bins.csv',np.column_stack([edges[:-1],edges[1:],count,np.sqrt(count)]),delimiter=',',header='mass_low_GeV,mass_high_GeV,count,sqrt_count',comments='')
identity=json.loads((repo/'plotER/results/pb18_an_score/metadata/run2_score_identity.json').read_text());assert identity['matched']==10000 and identity['score_exact_matches']['BDT']==10000
meta={'source':str(src),'source_size':st.st_size,'source_mtime_ns':st.st_mtime_ns,'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'reference':str(ref),'reference_sha256':hashlib.sha256(ref.read_bytes()).hexdigest(),'reference_pdf_page':16,'reference_printed_page':14,'reference_table':3,'selection':'Score > 0.06 && BQvalue < 0.13','score_correspondence':'10000 evenly spaced REflated candidates matched by exact Bmass/Bpt/By to original ntmix_flatten; Score equals BDT for all matched rows','input_entries':len(mask),'selected_entries':int(n),'histogram_entries':int(count.sum()),'underflow':int(h.GetBinContent(0)),'overflow':int(h.GetBinContent(81)),'mass_range':[3.6,4.0],'bin_width_GeV':.005,'weights':'unit DATA counts','scope':'Published BDT threshold plus user-authorized BQvalue < 0.13 on existing REflated DATA; not full AN selection reproduction','additional_an_cuts_not_applied':{'Bpt_range':[15,50],'abs_By_lt':1.6,'centrality_percent':[0,90],'abs_PVz_lt':15},'Q_definition_note':'REflated BQvalue = Bmass - 3.096916 - Btktkmass; AN text states m_mumupipi - m_mumu - m_pipi. Original Bmumumass is distinct from flat Bujmass.','verification':'Independent numpy vs ROOT selection count and every histogram bin agree; finite inputs; no fit or optimization','ROOT_version':R.gROOT.GetVersion(),'plot_template':'plotER/macros/plot_PbPb_Data_Bmass_local.C single-year layout, binning, color/fill and HIST style'}
assert src.stat().st_mtime_ns==st.st_mtime_ns;(out/'metadata/manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
(out/'README.md').write_text('''# PbPb18：AN2019 的最终 BDT cut

[质量谱 PDF](Bmass_score0p06.pdf)

AN2019_021_v10 的印刷页14（PDF第16页），第4.1节表3明确给出 **BDT > 0.06**。无需重新扫描 optimal。现有 REflated DATA 的 Score 已与原始旧 BDT 抽查10000条，全部逐值相同。

本图对现有 `legacy_run2/REflated/flat_ntmix_PbPb18_DATA.root:ntmix` 追加 `Score > 0.06 && BQvalue < 0.13`；输入共 %d 条，选后 **%d** 条。范围3.60–4.00 GeV，5 MeV分箱，DATA单位计数，复用已有单年份Bmass宏的布局与样式。没有拟合。

这是应用文献 score cut 后的现有样本质量谱。AN完整最终选择还包括15<Bpt<50、abs(By)<1.6、中心度0–90%%、abs(PVz)<15；本次没有追加这些条件，不宣称完整复现AN图20。

尤其注意Q定义：当前BQvalue用固定3.096916 GeV，AN文字写逐候选双muon质量。原始Bmumumass与flat Bujmass也不同，不能直接互相替换。原始flatten的hiBin与REflated CentBin前1000候选数值相同，不能未经核实自动除以2。本次依用户明确指令使用现有BQvalue<0.13；当前flat无Bmumumass或另一个已核验等价的Q变量，Bujmass不能直接替代原始Bmumumass。

[来源与计数验收](metadata/manifest.json) · [逐bin计数](metadata/Bmass_bins.csv) · [ROOT直方图](metadata/Bmass_score0p06.root) · [Score核验](metadata/run2_score_identity.json)
'''%(len(mask),n))
print(json.dumps(meta,indent=2))
