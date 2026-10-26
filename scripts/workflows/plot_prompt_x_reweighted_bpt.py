#!/usr/bin/env python3
from pathlib import Path
from array import array
import argparse,hashlib,json,math,shutil
import numpy as np
import uproot
import ROOT
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
assert not args.output.exists();args.output.mkdir(parents=True);metadata=args.output/'metadata';metadata.mkdir()
repo=Path(__file__).resolve().parents[2];root=repo/'X3872'
manifest=json.loads((root/'metadata/mc_manifest.json').read_text());assert manifest['status']=='PASS' and manifest['source5_policy']=='zero'
assert json.loads((root/'metadata/independent_acceptance.json').read_text())['status']=='PASS'
def digest(path):
 checksum=hashlib.sha256()
 with path.open('rb') as stream:
  for chunk in iter(lambda:stream.read(8388608),b''):checksum.update(chunk)
 return checksum.hexdigest()
samples={};provenance={};report={'selection':None,'weight':'pThatreweight','normalization':'unit area separately for each year','years':{}}
for record in manifest['files']:
 if record['particle']!='x3872' or record['kind']!='prompt':continue
 year=record['year'];path=Path(record['output']);path=path if path.is_absolute() else repo/path
 stat=path.stat();checksum=digest(path);assert checksum==record['output_sha256']
 with uproot.open(path) as file:values=file['ntmix_X3872'].arrays(['Bpt','pThatreweight','pThats','pthat'],library='np')
 weights=values['pThatreweight'];assert np.isfinite(weights).all() and (weights>=0).all();assert np.all(weights[(values['pThats']==5)|(values['pthat']<values['pThats'])]==0)
 samples[year]=values;provenance[year]={'path':str(path),'sha256':checksum,'size':stat.st_size,'mtime_ns':stat.st_mtime_ns}
 report['years'][year]={'entries':len(weights),'positive_entries':int((weights>0).sum()),'sumw':float(weights.sum()),'sumw2':float(weights@weights),'Bpt_range':[float(values['Bpt'].min()),float(values['Bpt'].max())]}
ROOT.gROOT.SetBatch(True);ROOT.gStyle.SetOptStat(0);ROOT.gStyle.SetTextFont(42)
maximum=5*math.ceil(max(values['Bpt'].max() for values in samples.values())/5);edges=np.arange(0,maximum+5,5,dtype=float)
canvas=ROOT.TCanvas('Bpt','',850,700);canvas.SetLeftMargin(.15);canvas.SetRightMargin(.04);canvas.SetBottomMargin(.13);canvas.SetTopMargin(.08);canvas.SetLogy()
rootfile=ROOT.TFile.Open(str(metadata/'histograms.root'),'RECREATE');histograms=[];histogram_report={}
for year,color in [('2018',ROOT.kBlue+1),('2023',ROOT.kOrange+7)]:
 values=samples[year];coordinates=values['Bpt'].astype(float);weights=values['pThatreweight'].astype(float);counts=np.histogram(coordinates,edges,weights=weights)[0];squares=np.histogram(coordinates,edges,weights=weights**2)[0]
 histogram=ROOT.TH1D('Bpt_'+year,'',len(edges)-1,array('d',edges));histogram.Sumw2();histogram.FillN(len(coordinates),np.ascontiguousarray(coordinates),np.ascontiguousarray(weights))
 assert histogram.GetBinContent(0)==histogram.GetBinContent(len(edges))==0
 assert np.allclose([histogram.GetBinContent(index+1) for index in range(len(counts))],counts,rtol=1e-8,atol=1e-7)
 assert np.allclose([histogram.GetBinError(index+1)**2 for index in range(len(counts))],squares,rtol=1e-8,atol=1e-7)
 histogram.Write();histogram.Scale(1/weights.sum());histogram.SetName('Bpt_'+year+'_unit_area');histogram.SetLineColor(color);histogram.SetMarkerColor(color);histogram.SetLineWidth(2);histogram.SetMarkerStyle(20 if year=='2018' else 24);histogram.SetMarkerSize(.65);histogram.Write();histograms.append(histogram)
 histogram_report[year]={'edges':edges.tolist(),'sumw':counts.tolist(),'sumw2':squares.tolist(),'normalization':float(weights.sum()),'no_additional_selection':True,'underflow':0,'overflow':0}
frame=histograms[0];frame.SetMinimum(1e-8);frame.SetMaximum(1.);frame.GetXaxis().SetTitle('RECO X(3872) p_{T} [GeV/c]');frame.GetYaxis().SetTitle('Unit-area weighted candidates / 5 GeV/c');frame.GetYaxis().SetTitleOffset(1.6)
frame.Draw('E1 HIST');histograms[1].Draw('E1 HIST SAME')
legend=ROOT.TLegend(.56,.66,.93,.83);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.033)
for histogram,year in zip(histograms,['2018','2023']):legend.AddEntry(histogram,year+' prompt X MC (weighted)','lp')
legend.Draw();latex=ROOT.TLatex();latex.SetNDC();latex.SetTextFont(42);latex.SetTextSize(.035);latex.DrawLatex(.15,.95,'PbPb prompt X(3872) MC');latex.SetTextSize(.029);latex.DrawLatex(.43,.60,'Existing flat RECO; no additional cuts');latex.DrawLatex(.43,.55,'Each year normalized to unit area');latex.DrawLatex(.43,.50,'Source 5 excluded; recomputed weights')
canvas.RedrawAxis();canvas.Print(str(args.output/'Bpt_weighted.pdf'));rootfile.Close()
report['histograms']=histogram_report;report['audit']={'ROOT_numpy_bins_and_sumw2':'PASS','input_manifest_hashes':'PASS','source5_and_below_threshold_zero':'PASS'}
for year,record in provenance.items():
 stat=Path(record['path']).stat();assert (stat.st_size,stat.st_mtime_ns)==(record['size'],record['mtime_ns'])
shutil.copy2(__file__,metadata/Path(__file__).name)
provenance['script_sha256']=digest(Path(__file__));provenance['ROOT_version']=ROOT.gROOT.GetVersion()
(metadata/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n');(metadata/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
(args.output/'README.md').write_text("# 新权重prompt X：2018与2023的Bpt分布\n\n[对比图](Bpt_weighted.pdf)。使用X3872下两年份新重算prompt X MC的全部RECO候选，不加Bpt、By、BQvalue、Bmass或ML筛选。来源5条目及低于自身来源阈值的条目权重为0，不贡献分布。\n\n沿用旧图分箱与样式：5 GeV/bin、0–255 GeV、对数纵轴；每年份分别归一化到单位面积，误差为sqrt(sum(w*w))/总权重。既有flatten选择不同（2018最低Bpt约14.56、2023约6.33 GeV），低端差异包含此影响。\n\n输入SHA与已验收新样本manifest一致；ROOT/NumPy逐bin加权计数和误差平方一致，无上下溢出。见[计数与验收](metadata/summary.json)、[输入来源](metadata/provenance.json)。原含来源5图保留在../prompt_x_pthat5/。\n")
print(json.dumps(report['years'],indent=2));print('PASS')
