#!/usr/bin/env python3
from pathlib import Path
from array import array
import argparse,hashlib,json,math,shutil
import numpy as np
import uproot
import ROOT
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
assert not args.output.exists();args.output.mkdir(parents=True);metadata=args.output/'metadata';metadata.mkdir()
repo=Path(__file__).resolve().parents[2]
previous=repo/'plotER/results/prompt_x_pthat5/metadata/provenance.json'
inputs=json.loads(previous.read_text())
def digest(path):
 checksum=hashlib.sha256()
 with path.open('rb') as stream:
  for chunk in iter(lambda:stream.read(8388608),b''):checksum.update(chunk)
 return checksum.hexdigest()
samples={};provenance={};report={'selection':'Bpt > 15 && Bpt < 50 && abs(By) < 1.6 && BQvalue < 0.15','weight':'pThatreweight','normalization':'unit area separately for each year','years':{}}
for year in ['2018','2023']:
 record=inputs[year];path=Path(record['path'])
 stat=path.stat();checksum=digest(path);assert checksum==record['sha256'], 'Original input changed since prompt_x_pthat5'
 with uproot.open(path) as file:values=file['ntmix_X3872'].arrays(['Bpt','By','BQvalue','pThatreweight'],library='np')
 weights=values['pThatreweight'];assert np.isfinite(weights).all() and (weights>=0).all()
 selected=(values['Bpt']>15)&(values['Bpt']<50)&(abs(values['By'])<1.6)&(values['BQvalue']<.15)
 check=ROOT.TFile.Open(str(path));assert check.Get('ntmix_X3872').GetEntries(report['selection'])==int(selected.sum());check.Close()
 values={key:val[selected] for key,val in values.items()};weights=values['pThatreweight']
 samples[year]=values;provenance[year]={'path':str(path),'sha256':checksum,'size':stat.st_size,'mtime_ns':stat.st_mtime_ns}
 report['years'][year]={'entries':len(weights),'positive_entries':int((weights>0).sum()),'sumw':float(weights.sum()),'sumw2':float(weights@weights),'Bpt_range':[float(values['Bpt'].min()),float(values['Bpt'].max())]}
ROOT.gROOT.SetBatch(True);ROOT.gStyle.SetOptStat(0);ROOT.gStyle.SetTextFont(42)
edges=np.arange(10,105,5,dtype=float)
canvas=ROOT.TCanvas('Bpt','',850,850)
top=ROOT.TPad('top','',0,.30,1,1);bottom=ROOT.TPad('bottom','',0,0,1,.30)
for pad in [top,bottom]:pad.SetLeftMargin(.15);pad.SetRightMargin(.04);pad.Draw()
top.SetTopMargin(.09);top.SetBottomMargin(.02);top.SetLogy();bottom.SetTopMargin(.025);bottom.SetBottomMargin(.32);top.cd()
rootfile=ROOT.TFile.Open(str(metadata/'histograms.root'),'RECREATE');histograms=[];histogram_report={}
for year,color in [('2018',ROOT.kBlue+1),('2023',ROOT.kOrange+7)]:
 values=samples[year];coordinates=values['Bpt'].astype(float);weights=values['pThatreweight'].astype(float);counts=np.histogram(coordinates,edges,weights=weights)[0];squares=np.histogram(coordinates,edges,weights=weights**2)[0]
 histogram=ROOT.TH1D('Bpt_'+year,'',len(edges)-1,array('d',edges));histogram.Sumw2();histogram.FillN(len(coordinates),np.ascontiguousarray(coordinates),np.ascontiguousarray(weights))
 assert histogram.GetBinContent(0)==histogram.GetBinContent(len(edges))==0
 assert np.allclose([histogram.GetBinContent(i+1) for i in range(len(counts))],counts,rtol=1e-8,atol=1e-7)
 assert np.allclose([histogram.GetBinError(i+1)**2 for i in range(len(counts))],squares,rtol=1e-8,atol=1e-7)
 histogram.Write();histogram.Scale(1/weights.sum());histogram.SetName('Bpt_'+year+'_unit_area');histogram.SetLineColor(color);histogram.SetMarkerColor(color);histogram.SetLineWidth(2);histogram.SetMarkerStyle(20 if year=='2018' else 24);histogram.SetMarkerSize(.7);histogram.Write();histograms.append(histogram)
 histogram_report[year]={'edges':edges.tolist(),'sumw':counts.tolist(),'sumw2':squares.tolist(),'normalization':float(weights.sum()),'underflow':0,'overflow':0}
frame=histograms[0];frame.SetMinimum(1e-3);frame.SetMaximum(1.);frame.GetXaxis().SetLabelSize(0);frame.GetYaxis().SetTitle('Unit-area weighted candidates / 5 GeV/c');frame.GetYaxis().SetTitleOffset(1.5)
frame.Draw('E1 HIST');histograms[1].Draw('E1 HIST SAME')
legend=ROOT.TLegend(.53,.70,.94,.86);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.032)
for histogram,year in zip(histograms,['2018','2023']):legend.AddEntry(histogram,year+' prompt X MC (weighted)','lp')
legend.Draw();latex=ROOT.TLatex();latex.SetNDC();latex.SetTextFont(42);latex.SetTextSize(.035);latex.DrawLatex(.15,.95,'PbPb prompt X(3872) MC');latex.SetTextSize(.031)
for y,label in [(.61,'15 < p_{T} < 50 GeV/c, |y| < 1.6'),(.55,'BQvalue < 0.15 GeV/c^{2}'),(.49,'Each year normalized to unit area'),(.43,'Source 5 included; original weights')]:latex.DrawLatex(.50,y,label)
top.RedrawAxis();bottom.cd()
ratio_frame=ROOT.TH1D('ratio_frame','',len(edges)-1,array('d',edges));ratio_frame.SetMinimum(0.);ratio_frame.SetMaximum(2.)
ratio_frame.GetXaxis().SetTitle('RECO X(3872) p_{T} [GeV/c]');ratio_frame.GetYaxis().SetTitle('2023 / 2018')
ratio_frame.GetXaxis().SetTitleSize(.12);ratio_frame.GetXaxis().SetLabelSize(.10);ratio_frame.GetXaxis().SetTitleOffset(1.05);ratio_frame.GetYaxis().SetTitleSize(.10);ratio_frame.GetYaxis().SetLabelSize(.09);ratio_frame.GetYaxis().SetTitleOffset(.58);ratio_frame.GetYaxis().SetNdivisions(505);ratio_frame.Draw('AXIS')
graph=ROOT.TGraphErrors();graph.SetName('ratio_2023_over_2018');graph.SetMarkerStyle(20);graph.SetMarkerSize(.8);ratio_records=[]
for i in range(1,len(edges)):
 a=histograms[1].GetBinContent(i);b=histograms[0].GetBinContent(i)
 if b<=0:continue
 ea=histograms[1].GetBinError(i);eb=histograms[0].GetBinError(i);ratio=a/b;error=math.sqrt((ea/b)**2+(a*eb/b**2)**2)
 n=graph.GetN();graph.SetPoint(n,float((edges[i-1]+edges[i])/2),ratio);graph.SetPointError(n,0,error)
 a_np=histogram_report['2023']['sumw'][i-1]/histogram_report['2023']['normalization'];b_np=histogram_report['2018']['sumw'][i-1]/histogram_report['2018']['normalization'];va=histogram_report['2023']['sumw2'][i-1]/histogram_report['2023']['normalization']**2;vb=histogram_report['2018']['sumw2'][i-1]/histogram_report['2018']['normalization']**2
 assert np.isclose(ratio,a_np/b_np,rtol=1e-8);assert np.isclose(error,math.sqrt(va/b_np**2+a_np**2*vb/b_np**4),rtol=1e-8)
 ratio_records.append({'lower':float(edges[i-1]),'upper':float(edges[i]),'ratio':ratio,'error':error})
assert len(ratio_records)==7
line=ROOT.TLine(10,1,100,1);line.SetLineStyle(2);line.SetLineColor(ROOT.kGray+2);line.Draw();graph.Draw('P SAME');graph.Write();bottom.RedrawAxis();canvas.cd();canvas.Print(str(args.output/'Bpt_weighted.pdf'));rootfile.Close()
report['histograms']=histogram_report;report['ratio']=ratio_records;report['ratio_error']='Independent-year propagation of sumw2 errors with normalization factors held fixed, matching main panel';report['audit']={'ROOT_numpy_bins_and_sumw2':'PASS','ROOT_selection_counts':'PASS','ratio_numpy_values_and_errors':'PASS','original_prompt_x_pthat5_input_hashes':'PASS'}
for year,record in provenance.items():
 stat=Path(record['path']).stat();assert (stat.st_size,stat.st_mtime_ns)==(record['size'],record['mtime_ns'])
shutil.copy2(__file__,metadata/Path(__file__).name);provenance['script_sha256']=digest(Path(__file__));provenance['ROOT_version']=ROOT.gROOT.GetVersion()
(metadata/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n');(metadata/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
(args.output/'README.md').write_text("# 原含来源5的prompt X：指定fiducial内Bpt与ratio\n\n[对比图](Bpt_weighted.pdf)。2018、2023均使用此前prompt_x_pthat5对应的Henrique原prompt X MC及现成pThatreweight，输入SHA与当时记录一致，选择15 < Bpt < 50、abs(By) < 1.6、BQvalue < 0.15，不加其它筛选。保留来源5及其原权重贡献，不使用本项目新重算权重。\n\n横轴10–100 GeV、5 GeV/bin，主图对数纵轴，每年份在选后分别归一化到单位面积。误差沿用sqrt(sum(w*w))/总权重。下图为2023/2018的归一化分布比值，误差按两年份独立传播主图误差（归一化因子固定，未计归一化引入的bin间相关）；分母为0的区间不画点。Bpt<50的筛选使50–100 GeV区域为空。\n\n输入SHA、ROOT独立选择计数、ROOT/NumPy逐bin及sumw2、ratio数值/误差核验通过。见[完整数值](metadata/summary.json)和[输入来源](metadata/provenance.json)。\n")
print(json.dumps({'years':report['years'],'ratio':ratio_records},indent=2));print('PASS')
