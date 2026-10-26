#!/usr/bin/env python3
# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
"""PbPb23 prompt/nonprompt signed BLxy comparison, common v27 precut, before ML."""
from pathlib import Path
import argparse,json,hashlib,shutil,subprocess
import numpy as np,uproot,ROOT as R
ap=argparse.ArgumentParser();ap.add_argument('--output',required=True,type=Path);args=ap.parse_args();out=args.output.resolve();out.mkdir(parents=True,exist_ok=False);art=out/'artifacts';art.mkdir()
assert R.gROOT.GetVersion()=='6.32.02';R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0)
base=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23')
selection='Bpt > 15 && Bpt < 50 && abs(By) < 1.6 && BQvalue < 0.15 && Btrk2dR > 0 && Btrk2dR < 0.35'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 return h.hexdigest()
meta={'selection':selection,'ml_cut':None,'mass_cut':None,'lxy_cut':None,'variable':'signed BLxy [cm], not abs or pseudo-proper length','weight':'pThatreweight','normalization':'full selected finite support separately per sample, including displayed underflow/overflow','marker_cm':.01,'root_version':R.gROOT.GetVersion(),'species':{}}
root=R.TFile(str(art/'histograms.root'),'RECREATE')
for species,label,slug in [('PSI2S','#psi(2S)','psi2s'),('X3872','X(3872)','x')]:
 samples={};info={}
 for kind,suffix in [('prompt',''),('nonprompt','_nonPrompt')]:
  p=base/f'flat_ntmix_PbPb23_MC_{species}{suffix}.root';digest=sha(p)
  with uproot.open(p) as f:a=f['ntmix_'+species].arrays(['BLxy','pThatreweight','Bpt','By','BQvalue','Btrk2dR'],library='np')
  mask=(a['Bpt']>15)&(a['Bpt']<50)&(abs(a['By'])<1.6)&(a['BQvalue']<.15)&(a['Btrk2dR']>0)&(a['Btrk2dR']<.35)
  assert np.all(np.isfinite(a['BLxy'][mask])) and np.all(np.isfinite(a['pThatreweight'][mask])) and np.all(a['pThatreweight'][mask]>0)
  x=a['BLxy'][mask].astype(float);w=a['pThatreweight'][mask].astype(float);samples[kind]=(x,w)
  info[kind]={'path':str(p),'sha256':digest,'source_entries':len(mask),'selected_entries':len(x),'sumw':float(w.sum()),'sumw2':float(w@w),'fraction_below_0p01':float(w[x<.01].sum()/w.sum()),'raw_below_0p01':int((x<.01).sum()),'views':{}}
  assert sha(p)==digest
 c=R.TCanvas('canvas_'+slug,'',1200,600);c.Divide(2,1);keep=[]
 for panel,(view,nb,lo,hi,logy) in enumerate([('near_zero',75,-.05,.10,False),('wide',130,-.10,2.50,True)],1):
  pad=c.cd(panel);pad.SetLeftMargin(.14);pad.SetRightMargin(.04);pad.SetBottomMargin(.13);pad.SetTopMargin(.10);pad.SetLogy(logy);hs=[]
  for kind,color,marker in [('prompt',R.kOrange+7,24),('nonprompt',R.kBlue+1,20)]:
   x,w=samples[kind];edges=np.linspace(lo,hi,nb+1);inside=(x>=lo)&(x<hi);bins=np.histogram(x[inside],edges,weights=w[inside])[0];v=np.histogram(x[inside],edges,weights=w[inside]**2)[0];h=R.TH1D(slug+'_'+kind+'_'+view,'',nb,lo,hi);h.SetDirectory(0);h.Sumw2()
   for i in range(nb):h.SetBinContent(i+1,bins[i]/w.sum());h.SetBinError(i+1,np.sqrt(v[i])/w.sum())
   # Independent ROOT selection and histogram check.
   f=R.TFile.Open(info[kind]['path']);t=f.Get('ntmix_'+species);R.gROOT.cd();check=R.TH1D('check_hist','',nb,lo,hi);check.Sumw2();t.Draw('BLxy>>check_hist','('+selection+')*pThatreweight','goff')
   assert t.GetEntries(selection)==len(x)
   for i in range(nb):assert abs(check.GetBinContent(i+1)/w.sum()-h.GetBinContent(i+1))<1e-10 and abs(check.GetBinError(i+1)/w.sum()-h.GetBinError(i+1))<1e-10
   check.SetDirectory(0);del check;f.Close();pad.cd()
   h.SetLineColor(color);h.SetMarkerColor(color);h.SetMarkerStyle(marker);h.SetMarkerSize(.35);h.SetLineWidth(2);hs.append(h)
   info[kind]['views'][view]={'bins':nb,'range':[lo,hi],'underflow_weight_fraction':float(w[x<lo].sum()/w.sum()),'overflow_weight_fraction':float(w[x>=hi].sum()/w.sum()),'bin_contents':(bins/w.sum()).tolist(),'bin_errors':(np.sqrt(v)/w.sum()).tolist()}
   root.cd();h.Write();pad.cd()
  high=max(h.GetBinContent(i)+h.GetBinError(i) for h in hs for i in range(1,nb+1));ymin=1e-6 if logy else 0.;ymax=high*(30 if logy else 1.55)
  hs[0].SetMinimum(ymin);hs[0].SetMaximum(ymax);hs[0].GetXaxis().SetTitle('L_{xy} [cm]');hs[0].GetYaxis().SetTitle('Weighted fraction / bin');hs[0].GetYaxis().SetTitleOffset(1.6);hs[0].Draw('HIST E');hs[1].Draw('HIST E SAME')
  line=R.TLine(.01,ymin,.01,ymax);line.SetLineStyle(2);line.SetLineWidth(2);line.SetLineColor(R.kBlack);line.Draw()
  leg=R.TLegend(.58,.68,.94,.88);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.031);leg.AddEntry(hs[0],'Prompt MC','le');leg.AddEntry(hs[1],'Non-prompt MC','le');leg.AddEntry(line,'L_{xy} = 0.01 cm','l');leg.Draw()
  text=R.TLatex();text.SetNDC();text.SetTextSize(.037);text.DrawLatex(.16,.94,'PbPb23 '+label);text.SetTextSize(.027);text.DrawLatex(.16,.86,'v27 fiducial, before ML cut');text.DrawLatex(.16,.81,'15 < p_{T} < 50 GeV/c, |y| < 1.6')
  keep.extend([*hs,line,leg,text])
 c.SaveAs(str(out/(slug+'_lxy.pdf')));c.Close();meta['species'][slug]=info
root.Close();subprocess.run(['pdfunite',str(out/'psi2s_lxy.pdf'),str(out/'x_lxy.pdf'),str(out/'validation_all.pdf')],check=True)
meta['audit']='PASS: all histogram bins and Sumw2 independently reproduced with ROOT TTree::Draw; selected counts and input hashes checked'
(out/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n');shutil.copy2(__file__,art/Path(__file__).name)
print(json.dumps({s:{k:r['fraction_below_0p01'] for k,r in d.items()} for s,d in meta['species'].items()},indent=2))
