from pathlib import Path
import ROOT as R,numpy as np,json,hashlib,shutil,os
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0)
o=Path('/eos/home-l/leyao/pbpb_work/X_analysis/Analysis_CODES/plotER/results/pb18_an_score');src=o/'metadata/Bmass_score0p06.root';sha=hashlib.sha256(src.read_bytes()).hexdigest();pdf=o/'Bmass_score0p06_10MeV.pdf';assert not pdf.exists()
f=R.TFile.Open(str(src));h5=f.Get('Bmass_score0p06');assert h5.GetNbinsX()==80;h=h5.Rebin(2,'Bmass_score0p06_10MeV');h.SetDirectory(0);f.Close()
a=np.genfromtxt(o/'metadata/Bmass_bins.csv',delimiter=',',names=True);expected=a['count'].reshape(40,2).sum(axis=1);assert np.array_equal(expected,[h.GetBinContent(i) for i in range(1,41)]);assert int(h.Integral())==19699
h.GetYaxis().SetTitle('Entries / 10 MeV/c^{2}')
c=R.TCanvas('cPb18AnScore10','',1000,650);c.SetLeftMargin(.10);c.SetRightMargin(.36);c.SetTopMargin(.06);c.SetBottomMargin(.12)
h.SetLineColor(R.kBlue+1);h.SetLineWidth(2);h.SetFillColor(R.kBlue);h.SetFillStyle(3358);h.SetMinimum(0);h.SetMaximum(1.25*max(h.GetBinContent(i) for i in range(1, h.GetNbinsX()+1)));h.Draw('HIST')
l=R.TLegend(.66,.55,.985,.91);l.SetBorderSize(0);l.SetFillColor(R.kWhite);l.SetTextSize(.032);l.SetHeader('#bf{PbPb 2018 Data}');l.AddEntry(h,'N = 19,699','f')
for label in ['Existing REflated sample','Score (BDT) > 0.06','AN-19-021, Table 3','BQvalue < 0.13 GeV/c^{2}']:l.AddEntry(R.nullptr,label,'')
l.Draw();c.SaveAs(str(pdf))
root=o/'metadata/Bmass_score0p06_10MeV.root';assert not root.exists();f=R.TFile.Open(str(root),'RECREATE');h.Write();f.Close()
np.savetxt(o/'metadata/Bmass_bins_10MeV.csv',np.column_stack([a['mass_low_GeV'][::2],a['mass_high_GeV'][1::2],expected,np.sqrt(expected)]),delimiter=',',header='mass_low_GeV,mass_high_GeV,count,sqrt_count',comments='')
meta={'source_histogram':str(src),'source_sha256':sha,'selection':'Score > 0.06 && BQvalue < 0.13','bin_width_GeV':.010,'bins':40,'range':[3.6,4.0],'entries':19699,'verification':'All 40 bins equal sums of adjacent original 5 MeV CSV bins; original histogram unchanged','pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()};assert hashlib.sha256(src.read_bytes()).hexdigest()==sha;(o/'metadata/rebin_10MeV.json').write_text(json.dumps(meta,indent=2)+'\n');shutil.copy2('/tmp/pb18_10mev.py',o/'metadata/pb18_10mev.py')
p=o/'README.md';old=p.read_text();anchor='[质量谱 PDF](Bmass_score0p06.pdf)';assert old.count(anchor)==1;new=old.replace(anchor,'[5 MeV/bin 图](Bmass_score0p06.pdf) · [10 MeV/bin 图](Bmass_score0p06_10MeV.pdf)',1);new+='\n10 MeV/bin 图由已验收的5 MeV直方图相邻两bin相加，选后总数仍为19699，逐bin一致性通过。[合并记录](metadata/rebin_10MeV.json)。\n';backup=o/'metadata/originals/README_before_10MeV.md';assert not backup.exists();shutil.copy2(p,backup);tmp=p.with_name('.README.10MeV.tmp');tmp.write_text(new);shutil.copymode(p,tmp);assert p.read_text()==old;os.replace(tmp,p)
print('PASS: 40 bins, 19699 candidates; original 5 MeV histogram retained')
