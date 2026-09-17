from pathlib import Path
import json,csv,re,hashlib,subprocess,argparse
import numpy as np,uproot,ROOT as R
ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--variable',choices=['Btktkmass']);args=ap.parse_args();repo=Path.cwd();out=Path(args.output);out.mkdir(parents=True,exist_ok=False);(out/'figures').mkdir();(out/'artifacts').mkdir()
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0)
src=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23');paths={'data':src/'flat_ntmix_PbPb23_DATA.root','mc':src/'flat_ntmix_PbPb23_MC_X3872.root'}
def fingerprint(p):
 st=p.stat();return {'path':str(p),'bytes':st.st_size,'mtime_ns':st.st_mtime_ns}
provenance={k:fingerprint(p) for k,p in paths.items()}
base=repo/'plotER/Validation/results/ppRef_X_r5_splot_ppref_snapshot_v1'
names=[n for n in json.loads((base/'artifacts/ppref_x_signed_sweight_all_common.json').read_text())['variables'] if n not in ['Bmass','signal_sWeight']];assert len(names)==67
if args.variable:names=[args.variable]
windows=[[3.75,3.85],[3.9,4.0]]
samples={};weights={};cutflow={}
for key,tree in [('mc','ntmix_X3872'),('data','ntmix')]:
 chunks={n:[] for n in names};wc=[];total=0;precount=0;passed=0
 with uproot.open(paths[key]) as f:
  tr=f[tree];assert set(names).issubset(tr.keys());read=sorted(set(names+['Bmass','Bpt','By','BQvalue']+(['pThatreweight'] if key=='mc' else [])))
  for a in tr.iterate(read,step_size=100000,library='np'):
   pre=(a['Bpt']>10)&(a['Bpt']<50)&(abs(a['By'])<1.6)&(a['BQvalue']<.15);total+=len(pre);precount+=int(pre.sum());mask=pre.copy()
   if key=='data':mask &= ((a['Bmass']>3.75)&(a['Bmass']<3.85))|((a['Bmass']>3.9)&(a['Bmass']<4.0))
   passed+=int(mask.sum())
   for n in names:chunks[n].append(a[n][mask])
   w=a['pThatreweight'][mask].astype(float) if key=='mc' else np.ones(int(mask.sum()));assert np.isfinite(w).all() and (w>0).all();wc.append(w)
 label='sideband' if key=='data' else 'mc';samples[label]={n:np.concatenate(v) for n,v in chunks.items()};weights[label]=np.concatenate(wc);cutflow[label]={'total':total,'preselection':precount,'selected':passed};print('SELECTED',label,cutflow[label],flush=True)
old={}
for n in names:
 vals=np.r_[samples['mc'][n],samples['sideband'][n]];vals=vals[np.isfinite(vals)];old[n]={'support_min':float(vals.min()),'support_max':float(vals.max())}
config={}
for line in (repo/'plotER/Validation/aux.h').read_text().splitlines():
 q=re.search(r'\{"([^"]+)",\s*"([^"]+)",\s*kNBins,\s*([^,]+),\s*([^,]+),\s*(true|false)',line)
 if q:
  try:config[q[1]]=(float(q[3]),float(q[4]),q[5]=='true',q[2].split(';')[1])
  except ValueError:pass
R.gInterpreter.Declare('#include "'+str(repo/'plotER/Validation/aux.h')+'"')
for q in R.getSignalVars('ntmix_X3872'):config[str(q.expr)]=(float(q.xmin),float(q.xmax),bool(q.absVal),str(q.title).split(';')[1])
if 'Btktkmass' in names:
 lo=np.floor(old['Btktkmass']['support_min']/.05)*.05;hi=(np.floor(old['Btktkmass']['support_max']/.05)+1)*.05
 config['Btktkmass']=(float(lo),float(hi),False,'m_{#pi#pi} [GeV/c^{2}]')
def cdf(x,w,y,v):
 z=np.union1d(x,y);ix=np.argsort(x,kind='stable');iy=np.argsort(y,kind='stable');a=np.r_[0,np.cumsum(w[ix])][np.searchsorted(x[ix],z,side='right')]/w.sum();b=np.r_[0,np.cumsum(v[iy])][np.searchsorted(y[iy],z,side='right')]/v.sum();return float(np.max(abs(a-b)))
metrics=[];root=R.TFile(str(out/'artifacts/histograms.root'),'RECREATE')
for j,n in enumerate(names):
 lo,hi,ab,title=config.get(n,(float(old[n]['support_min']),float(old[n]['support_max']),False,n))
 if hi<=lo:lo,hi=lo-.5,hi+.5
 edges=np.linspace(lo,hi,16);hs={};vals={};rec={'variable':n,'range':[lo,hi],'absolute_display':ab,'samples':{}}
 for k in samples:
  x=samples[k][n].astype(float);w=weights[k];ok=np.isfinite(x)&np.isfinite(w);x=x[ok];w=w[ok];vals[k]=(x,w);xx=abs(x) if ab else x;inside=(xx>=lo)&(xx<hi);hval=np.histogram(xx[inside],edges,weights=w[inside])[0];err=np.sqrt(np.histogram(xx[inside],edges,weights=w[inside]**2)[0]);norm=hval.sum();fallback=norm<=0
  if fallback:norm=w.sum()
  assert norm>0
  h=R.TH1D(n+'_'+k,'',15,lo,hi);h.SetDirectory(0);h.Sumw2()
  for i in range(15):h.SetBinContent(i+1,hval[i]/norm);h.SetBinError(i+1,err[i]/norm)
  hs[k]=h;rec['samples'][k]={'entries':len(x),'nonfinite':int((~ok).sum()),'sumw':float(w.sum()),'inrange_weight_fraction':float(hval.sum()/w.sum()),'normalization':'full support fallback' if fallback else 'inrange'}
 rec['cdf_sideband_mc']=cdf(*vals['sideband'],*vals['mc'])
 spH,mcH,sbH=hs['sideband'],hs['mc'],hs['sideband']
 for h,color,marker in [(sbH,R.kBlue+1,20),(mcH,R.kOrange+7,1)]:h.SetLineColor(color);h.SetMarkerColor(color);h.SetMarkerStyle(marker);h.SetLineWidth(2)
 band=mcH.Clone(n+'_band');band.SetDirectory(0);band.SetFillColorAlpha(R.kOrange+7,.3);band.SetMarkerSize(0)
 c=R.TCanvas('c'+str(j),'',760,650);top=R.TPad('top'+str(j),'',0,.30,1,1);bot=R.TPad('bot'+str(j),'',0,0,1,.30)
 for pad in [top,bot]:pad.SetLeftMargin(.14);pad.SetRightMargin(.04);pad.Draw()
 top.SetBottomMargin(.01);bot.SetTopMargin(.01);bot.SetBottomMargin(.33);top.cd();high=max(h.GetBinContent(i)+h.GetBinError(i) for h in hs.values() for i in range(1,16));low=min(0,min(h.GetBinContent(i)-h.GetBinError(i) for h in hs.values() for i in range(1,16)));spH.SetMaximum(max(high,.1)*1.7);spH.SetMinimum(low*1.12);spH.GetXaxis().SetLabelSize(0);spH.GetYaxis().SetTitle('Normalized entries');spH.Draw('E');band.Draw('E2 SAME');mcH.Draw('HIST SAME');sbH.Draw('E SAME');spH.Draw('E SAME')
 leg=R.TLegend(.57,.76,.94,.91);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.034);leg.AddEntry(mcH,'Prompt X MC','l');leg.AddEntry(sbH,'Sideband DATA (raw)','lep');leg.Draw();lab=R.TLatex();lab.SetNDC();lab.SetTextSize(.047);lab.DrawLatex(.18,.87,'PbPb23 X(3872)');lab.SetTextSize(.031);lab.DrawLatex(.18,.80,'Preselection, no ML cut');lab.DrawLatex(.18,.74,'SB: (3.75,3.85) + (3.90,4.00)');lab.DrawLatex(.18,.68,'CDF D(SB,MC) = %.3f'%rec['cdf_sideband_mc'])
 bot.cd();axis=R.TH1D(n+'_ratioaxis','',15,lo,hi);axis.SetDirectory(0);axis.SetMinimum(0);axis.SetMaximum(5);axis.GetYaxis().SetTitle('Sideband / MC');axis.GetYaxis().SetTitleSize(.09);axis.GetYaxis().SetLabelSize(.08);axis.GetYaxis().SetTitleOffset(.7);axis.GetYaxis().SetNdivisions(304);axis.GetXaxis().SetTitle(title);axis.GetXaxis().SetTitleSize(.11);axis.GetXaxis().SetLabelSize(.10);axis.Draw('AXIS');gs=[];clips=0
 for k,color,marker in [('sideband',R.kBlue+1,20)]:
  h=hs[k];g=R.TGraphErrors();g.SetName(n+'_'+k+'_ratio');g.SetLineColor(color);g.SetMarkerColor(color);g.SetMarkerStyle(marker)
  for i in range(1,16):
   den=mcH.GetBinContent(i)
   if den<=0:continue
   val=h.GetBinContent(i)/den;err=np.hypot(h.GetBinError(i)/den,h.GetBinContent(i)*mcH.GetBinError(i)/den**2);idx=g.GetN();g.SetPoint(idx,h.GetBinCenter(i),val);g.SetPointError(idx,0,float(err));clips+=int(val<0 or val>5)
  g.Draw('P SAME');gs.append(g)
 line=R.TLine(lo,1,hi,1);line.SetLineStyle(2);line.Draw();rec['ratio_points_outside_display']=clips
 if clips:lab.SetTextSize(.075);lab.DrawLatex(.18,.85,f'{clips} ratio points outside [0,5]')
 c.SaveAs(str(out/'figures'/f'{n}.pdf'));root.cd()
 for obj in [*hs.values(),*gs]:obj.Write()
 c.Close();metrics.append(rec)
root.Close();subprocess.run(['pdfunite',*[str(out/'figures'/f'{n}.pdf') for n in names],str(out/'figures/validation_all.pdf')],check=True)
for p in (out/'figures').glob('*.pdf'):assert 'Pages:' in subprocess.check_output(['pdfinfo',str(p)],text=True)
for k,p in paths.items():assert fingerprint(p)==provenance[k]
summary={'status':'PASS','selection':'Bpt > 10 && Bpt < 50 && abs(By) < 1.6 && BQvalue < 0.15','sideband_definition':'raw DATA, strict open mass intervals','sideband_windows':windows,'mc_mass_cut':None,'counts':{k:len(w) for k,w in weights.items()},'cutflow':cutflow,'sumw':{k:float(w.sum()) for k,w in weights.items()},'mc_weight_unique':np.unique(weights['mc']).tolist(),'variables':len(names),'input_provenance':provenance,'notes':['Prompt X MC pThatreweight; raw sideband DATA unit weight.','No ML, dR, or displacement significance cut; no fitting or subtraction.','Histogram normalized in visible range; full-support fallback for empty visible range recorded per variable.','CDF full finite raw support; separation diagnostic, not MC-to-signal validation.','Input file size/mtime checked before/after; no claim of full source SHA256.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(out/'variable_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n');print(json.dumps(summary,indent=2))
