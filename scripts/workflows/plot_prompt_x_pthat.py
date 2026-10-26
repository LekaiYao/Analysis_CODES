#!/usr/bin/env python3
"""Prompt X GEN pthat/Gpt comparison, AN2019_021 Figures 2/3 style.
2018 cache is built by the archived extraction script; verify it against both
original event branches and the complete existing GEN flat before plotting.
Run in LCG_106 (ROOT 6.32.02). Refuses to overwrite results.
"""
from pathlib import Path
import argparse,json,hashlib,shutil
from array import array
import numpy as np
import uproot
import ROOT as R

def dump(p,d): p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for s in iter(lambda:f.read(8388608),b''): h.update(s)
 return h.hexdigest()
def arrsha(a): return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def stat(p):
 s=p.stat();return {'path':str(p),'size':s.st_size,'mtime_ns':s.st_mtime_ns}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run2-cache',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
 out=args.output;assert not out.exists();out.mkdir(parents=True);(out/'pdf').mkdir();meta=out/'metadata';meta.mkdir()
 R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetTextFont(42)
 root=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872')
 paths={'run2_original':root/'legacy_run2/inputs/mc_prompt/prompt_x3872.root','run2_flat':root/'legacy_run2/REflated/flat_ntmix_PbPb18_MC_X3872.root','run3_flat':root/'PbPb23/flat_ntmix_PbPb23_MC_X3872.root'}
 provenance={k:stat(p) for k,p in paths.items()}
 with np.load(args.run2_cache) as f: a={k:f[k] for k in f.files}
 with uproot.open(paths['run2_original']) as f:
  ht=f['hiEvtAnalyzer/HiTree'];events=ht.arrays(['pthat','pthatweight','sample','run','evt'],library='np');ids=f['Bfinder/ntmix'].arrays(['RunNo','EvtNo'],library='np')
  assert np.array_equal(events['run'],ids['RunNo']) and np.array_equal(events['evt'],ids['EvtNo'])
  assert ht.num_entries==f['Bfinder/ntGen'].num_entries
  provenance['run2_original'].update(root_uuid=str(f.file.uuid),trees={'event':ht.object_path,'gen':f['Bfinder/ntGen'].object_path},event_entries=ht.num_entries,event_columns_sha256={k:arrsha(v) for k,v in events.items()},full_file_sha256=None)
 for k in ['pthat','sample','pthatweight']: assert np.array_equal(a[k],events[k][a['source_entry']])
 with uproot.open(paths['run2_flat']) as f:
  flat=f['ntGen'].arrays(['Gpt','Gy','GpdgId'],library='np')
  assert all(np.array_equal(a[k],flat[k]) for k in flat),'GEN flat must match the recovered signal candidates exactly'
  provenance['run2_flat'].update(root_uuid=str(f.file.uuid),tree='ntGen',sha256=sha(paths['run2_flat']))
 with uproot.open(paths['run3_flat']) as f:
  b=f['ntGen'].arrays(['Gpt','Gy','pthat','pThats','pThatreweight'],library='np');provenance['run3_flat'].update(root_uuid=str(f.file.uuid),tree='ntGen',sha256=sha(paths['run3_flat']))
 data={'2018':{'Gpt':a['Gpt'],'pthat':a['pthat'],'source':a['sample'].astype(np.int32),'weight':a['pthatweight'].astype(float)},'2023':{'Gpt':b['Gpt'],'pthat':b['pthat'],'source':b['pThats'],'weight':b['pThatreweight']}}
 assert (a['Gpt']>=5).all() and (abs(a['Gy'])<=2.4).all();assert (b['Gpt']>=5).all() and (abs(b['Gy'])<=2.4).all()
 shutil.copy2(args.run2_cache,meta/'pb18_prompt_x_gen.npz');shutil.copy2(__file__,meta/Path(__file__).name)
 provenance['run2_cache']={'sha256':sha(meta/'pb18_prompt_x_gen.npz'),'columns_sha256':{k:arrsha(v) for k,v in a.items()}}
 summary={'normalization':'Each curve separately normalized to sum of visible bins (0 <= x < 80); fraction per 2 GeV bin, not absolute yields. Overflow excluded from normalization and recorded.','selection':'Existing ntGen candidates: GisSignal==7, Gpt>=5, abs(Gy)<=2.4; 2018 recovered candidates exactly match existing flat GEN. No RECO, ML or additional cut.','pthat_counting':'One entry per selected GEN X candidate, including repeated event pthat if multiple selected X exist.','weights':{'2018':'hiEvtAnalyzer/HiTree.pthatweight, broadcast using ntGen event entry; no Ncoll or other weights','2023':'ntGen.pThatreweight, unchanged'},'source_labels':{'2018':'hiEvtAnalyzer/HiTree.sample','2023':'ntGen.pThats'},'years':{},'histograms':{}}
 edges=np.arange(0.,82.,2.);colors=[R.kAzure+1,R.kOrange+7,R.kGreen+2,R.kMagenta+1,R.kCyan+2];keep=[]
 rf=R.TFile.Open(str(meta/'histograms.root'),'RECREATE');book=out/'all_mc.pdf';canvas=R.TCanvas('c','',750,750);canvas.SetLeftMargin(.15);canvas.SetRightMargin(.04);canvas.SetBottomMargin(.13);canvas.SetTopMargin(.07);canvas.Print(str(book)+'[')
 for year,d in data.items():
  assert set(np.unique(d['source']))=={5,10,15,30,50};assert all(np.isfinite(v).all() for v in d.values());assert (d['weight']>=0).all()
  summary['years'][year]={'gen_candidates':len(d['Gpt']),'zero_weight_candidates':int((d['weight']==0).sum()),'sources':{str(s):int((d['source']==s).sum()) for s in [5,10,15,30,50]},'columns_sha256':{k:arrsha(v) for k,v in d.items()}}
  for variable in ['pthat','Gpt']:
   x=d[variable].astype(float);curves=[];stats={};name=f'pb{year}_{variable}'
   sets=[(f'pthat{s}',f'pThat{s} MC (unweighted)',d['source']==s,np.ones(len(x)),colors[i],1,2) for i,s in enumerate([5,10,15,30,50])]
   sets += [('total_unweighted','Total MC (unweighted)',np.ones(len(x),bool),np.ones(len(x)),R.kGray+2,2,2),('total_weighted','Total MC (weighted)',np.ones(len(x),bool),d['weight'],R.kRed+1,1,3)]
   for key,label,sel,w,color,style,width in sets:
    xx=x[sel];ww=w[sel];visible=(xx>=0)&(xx<80)
    counts=np.histogram(xx[visible],edges,weights=ww[visible])[0];sumw2=np.histogram(xx[visible],edges,weights=ww[visible]**2)[0];total=float(counts.sum());assert total>0
    # Independent ROOT binning and sum-of-squared-weights check.
    rh=R.TH1D(name+'_'+key+'_raw','',40,array('d',edges));rh.Sumw2();rh.FillN(len(xx),np.ascontiguousarray(xx,dtype='d'),np.ascontiguousarray(ww,dtype='d'))
    assert np.allclose([rh.GetBinContent(i+1) for i in range(40)],counts,rtol=1e-8,atol=1e-7)
    assert np.allclose([rh.GetBinError(i+1)**2 for i in range(40)],sumw2,rtol=1e-8,atol=1e-7)
    rh.Write();h=rh.Clone(name+'_'+key);h.Scale(1/total);h.SetLineColor(color);h.SetLineStyle(style);h.SetLineWidth(width);h.SetTitle('');h.Write();keep.extend([rh,h]);curves.append(h)
    stats[key]={'label':label,'entries':len(xx),'sumw_all':float(ww.sum()),'visible_sumw':total,'underflow_sumw':float(ww[xx<0].sum()),'overflow_sumw':float(ww[xx>=80].sum()),'raw_sumw':counts.tolist(),'raw_sumw2':sumw2.tolist(),'normalized':(counts/total).tolist()}
   assert np.array_equal(sum(np.array(stats['pthat'+str(s)]['raw_sumw']) for s in [5,10,15,30,50]),stats['total_unweighted']['raw_sumw'])
   canvas.Clear();canvas.SetLogy(True);frame=curves[-1];frame.SetMinimum(1e-7);frame.SetMaximum(1.)
   frame.GetXaxis().SetTitle('#hat{p}_{T} [GeV/c]' if variable=='pthat' else 'GEN X(3872) p_{T} [GeV/c]');frame.GetYaxis().SetTitle('Normalized candidates / 2 GeV/c')
   for axis in [frame.GetXaxis(),frame.GetYaxis()]:axis.SetTitleFont(42);axis.SetLabelFont(42);axis.SetTitleSize(.043);axis.SetLabelSize(.036)
   frame.GetYaxis().SetTitleOffset(1.65);frame.Draw('HIST')
   for h in curves[:-1]:h.Draw('HIST SAME')
   frame.Draw('HIST SAME')
   leg=R.TLegend(.57,.56,.94,.85);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextFont(42);leg.SetTextSize(.029)
   for h,st in zip(curves,stats.values()):leg.AddEntry(h,st['label'],'l')
   leg.Draw();keep.append(leg)
   tex=R.TLatex();tex.SetNDC();tex.SetTextFont(42);tex.SetTextSize(.036);tex.DrawLatex(.15,.947,f'PbPb {year}  |  Prompt X(3872) MC');tex.SetTextSize(.026);tex.DrawLatex(.18,.885,'GEN p_{T} #geq 5 GeV/c, |y| #leq 2.4');tex.DrawLatex(.18,.845,'Each curve: unit area in 0-80 GeV/c');keep.append(tex)
   canvas.RedrawAxis();canvas.Print(str(out/'pdf'/f'{name}.pdf'));canvas.Print(str(book));summary['histograms'][name]={'edges':edges.tolist(),'curves':stats}
 canvas.Print(str(book)+']');rf.Close()
 for k,p in paths.items():assert stat(p)=={n:provenance[k][n] for n in ['path','size','mtime_ns']},'Source changed'
 provenance['script_sha256']=sha(Path(__file__));provenance['root_version']=R.gROOT.GetVersion();summary['audit']={'status':'PASS','run2_event_ids_match':True,'run2_GEN_matches_flat_exactly':True,'numpy_ROOT_all_28_histograms_and_sumw2_match':True,'component_sum_equals_unweighted_total':True}
 dump(meta/'provenance.json',provenance);dump(meta/'distributions.json',summary);print(json.dumps({'audit':summary['audit'],'years':summary['years']},indent=2),flush=True)
if __name__=='__main__':main()
