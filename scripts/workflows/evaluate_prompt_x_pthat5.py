#!/usr/bin/env python3
from pathlib import Path
from array import array
import argparse,hashlib,json,math,shutil
import numpy as np
import uproot
import ROOT as ROOT

parser=argparse.ArgumentParser()
parser.add_argument('--run2-map',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
assert not args.output.exists()
args.output.mkdir(parents=True)
metadata=args.output/'metadata';metadata.mkdir()
root=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872')
repo=Path(__file__).resolve().parents[2]
paths={'2018':root/'legacy_run2/REflated/flat_ntmix_PbPb18_MC_X3872.root','2023':root/'PbPb23/flat_ntmix_PbPb23_MC_X3872.root'}

def digest(path):
 checksum=hashlib.sha256()
 with path.open('rb') as stream:
  for chunk in iter(lambda:stream.read(8388608),b''):checksum.update(chunk)
 return checksum.hexdigest()

def statistics(weights,denominator,events=None):
 weights=np.asarray(weights,dtype=float);total=float(weights.sum());square=float(weights@weights)
 result={'candidates':len(weights),'positive_candidates':int((weights>0).sum()),'sumw':total,'sumw2':square,'fraction_of_total_signal_weight':total/denominator if denominator else None,'candidate_neff':total*total/square if square else 0.,'largest_candidate_fraction_of_total':float(weights.max())/denominator if len(weights) and denominator else 0.,'top10_fraction_of_total':float(np.sort(weights)[-10:].sum())/denominator if denominator else None}
 if events is not None:
  unique,inverse=np.unique(events,return_inverse=True);sums=np.bincount(inverse,weights=weights);result.update(recorded_events=len(unique),positive_recorded_events=int((sums>0).sum()),event_neff=float(sums.sum()**2/(sums@sums)) if len(sums) and (sums@sums)>0 else 0.)
 else:result['recorded_events']=None
 return result

def selection(values):return (values['Bpt']>15)&(values['Bpt']<50)&(abs(values['By'])<1.6)&(values['BQvalue']<.15)

samples={};provenance={};report={'selection':'15 < Bpt < 50 && abs(By) < 1.6 && BQvalue < 0.15','MC_mass_cut':None,'background_sidebands':'3.80 < Bmass < 3.85 || 3.89 < Bmass < 3.94','years':{}}
for year,path in paths.items():
 before=path.stat();provenance[year]={'path':str(path),'size':before.st_size,'mtime_ns':before.st_mtime_ns,'sha256':digest(path)}
 with uproot.open(path) as file:
  columns=['Bpt','By','BQvalue','pThatreweight']+(['pThats','pthat'] if year=='2023' else [])
  values=file['ntmix_X3872'].arrays(columns,library='np')
 if year=='2018':
  with np.load(args.run2_map) as mapping:
   assert len(mapping['sample'])==len(values['Bpt'])
   values.update(pThats=mapping['sample'].copy(),pthat=mapping['pthat'].copy(),source_event=mapping['source_event'].copy())
 assert all(np.isfinite(value).all() for value in values.values()) and (values['pThatreweight']>=0).all()
 samples[year]=values;weights=values['pThatreweight'].astype(float);selected=selection(values);denominator=float(weights[selected].sum())
 def stat(mask):return statistics(weights[mask],denominator,values['source_event'][mask] if year=='2018' else None)
 masks={'all':selected,'source5':selected&(values['pThats']==5),'actual_pthat_le10':selected&(values['pthat']<=10),'source5_actual_pthat_le10':selected&(values['pThats']==5)&(values['pthat']<=10)}
 result={'flat_candidates':len(selected),'flat_Bpt_range':[float(values['Bpt'].min()),float(values['Bpt'].max())],'groups':{name:stat(mask) for name,mask in masks.items()},'by_source':{},'by_actual_pthat':{},'by_Bpt':{}}
 for source in [5,10,15,30,50]:result['by_source'][str(source)]=stat(selected&(values['pThats']==source))
 for lower,upper in [(0,10),(10,15),(15,30),(30,50),(50,np.inf)]:
  mask=selected&(values['pthat']>lower)&(values['pthat']<=upper);result['by_actual_pthat'][f'{lower}-{upper}']={'all':stat(mask),'source5':stat(mask&(values['pThats']==5))}
 for lower,upper in [(15,20),(20,25),(25,30),(30,40),(40,50)]:
  mask=selected&(values['Bpt']>=lower)&(values['Bpt']<upper);result['by_Bpt'][f'{lower}-{upper}']={'all':stat(mask),'source5':stat(mask&(values['pThats']==5))}
 assert sum(row['candidates'] for row in result['by_source'].values())==result['groups']['all']['candidates']
 assert np.isclose(sum(row['sumw'] for row in result['by_source'].values()),denominator)
 report['years'][year]=result

training=repo.parent/'XGBoost/output/training/fs_2v2_add1_X_pb23_v5_fid1_2v2_rw0_xgb_v1/X_pb23_v5_fid1_3v7_rw0_xgb_v1'
model_directory=repo.parent/'XGBoost/output/models'/training.parent.name/training.name
run_metadata=json.loads((model_directory/'run_metadata.json').read_text())
assert run_metadata['signal_selection']=='(Bpt > 15) and (Bpt < 50) and (abs(By) < 1.6) and (BQvalue < 0.15)'
assert run_metadata['background_selection']=='(((Bmass > 3.80) and (Bmass < 3.85)) or ((Bmass > 3.89) and (Bmass < 3.94))) and ((Bpt > 15) and (Bpt < 50) and (abs(By) < 1.6) and (BQvalue < 0.15))'
shutil.copy2(model_directory/'run_metadata.json',metadata/'ML_run_metadata.json')
shutil.copy2(model_directory/'model_config.json',metadata/'ML_model_config.json')
identity=json.loads((training/'split_identity.json').read_text());assert digest(training/'split_identity.npz')==identity['npz_sha256'];assert identity['sources'][0]['sha256']==provenance['2023']['sha256']
with np.load(training/'split_identity.npz') as split:
 assert np.array_equal(split['signal_entry'],np.flatnonzero(selection(samples['2023'])))
 signal_entries=split['signal_entry'];count=len(signal_entries)
 assert np.array_equal(np.sort(np.r_[split['train_index'],split['test_index']]),np.arange(count+len(split['background_entry'])))
 split_report={'tag':training.name,'background_candidates':len(split['background_entry']),'identity_sha256':identity['npz_sha256'],'groups':{}}
 for part in ['train','test']:
  chosen=split[part+'_index'];source_rows=signal_entries[chosen[chosen<count]];weights=samples['2023']['pThatreweight'][source_rows];source=samples['2023']['pThats'][source_rows];actual=samples['2023']['pthat'][source_rows];denominator=weights.sum()
  split_report['groups'][part]={name:statistics(weights[mask],denominator) for name,mask in {'all':np.ones(len(weights),bool),'source5':source==5,'source5_actual_pthat_le10':(source==5)&(actual<=10),'actual_pthat_le10':actual<=10}.items()}
roc=json.loads((training/'roc.json').read_text());split_report['existing_auc']={part:roc[part]['auc'] for part in ['train','test']};report['existing_training']=split_report
shutil.copy2(training/'split_identity.json',metadata/'split_identity.json');shutil.copy2(training/'split_identity.npz',metadata/'split_identity.npz')

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
legend.Draw();latex=ROOT.TLatex();latex.SetNDC();latex.SetTextFont(42);latex.SetTextSize(.035);latex.DrawLatex(.15,.95,'PbPb prompt X(3872) MC');latex.SetTextSize(.029);latex.DrawLatex(.43,.60,'Existing flat RECO; no additional cuts');latex.DrawLatex(.43,.55,'Each year normalized to unit area')
canvas.RedrawAxis();canvas.Print(str(args.output/'Bpt_weighted.pdf'));rootfile.Close()
report['histograms']=histogram_report;report['audit']={'ROOT_numpy_bin_counts_and_sumw2':'PASS','source_partition':'PASS','ML_actual_signal_selection_match':'PASS','independent_2023_event_identity':'UNAVAILABLE: flat lacks run/event identifiers','2018_event_mapping':'unique Bmass/Bpt/By match to original recorded event'}
for year,path in paths.items():
 after=path.stat();assert (after.st_size,after.st_mtime_ns)==(provenance[year]['size'],provenance[year]['mtime_ns'])
shutil.copy2(args.run2_map,metadata/'pb18_reco_source_map.npz');shutil.copy2(__file__,metadata/Path(__file__).name)
provenance['mapping_sha256']=digest(args.run2_map);provenance['script_sha256']=digest(Path(__file__));provenance['ROOT_version']=ROOT.gROOT.GetVersion()
(metadata/'provenance.json').write_text(json.dumps(provenance,indent=2)+'\n');(metadata/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
print(json.dumps({'years':{year:result['groups'] for year,result in report['years'].items()},'existing_training':split_report},indent=2),flush=True)
