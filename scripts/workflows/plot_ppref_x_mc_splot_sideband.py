from pathlib import Path
import json,csv,re,hashlib,subprocess,argparse
import numpy as np,uproot,ROOT as R
ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args();repo=Path.cwd();out=Path(args.output);out.mkdir(parents=True,exist_ok=False);(out/'figures').mkdir();(out/'artifacts').mkdir()
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0)
base=repo/'plotER/Validation/results/ppRef_X_r5_splot_ppref_snapshot_v1';src=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/ppRef24');paths={'data':src/'flat_ntmix_ppRef_DATA.root','mc':src/'flat_ntmix_ppRef_MC_X3872.root','splot':base/'artifacts/ppref_x_signed_sweight_all_common.root','model':repo/'fitER/ROOTfiles/ppRef_X_r5_fiducial_feasibility_ppref_snapshot_v1/nominalFitModel_ntmix_X3872_ppRef_X_r5_fiducial_feasibility_ppref_snapshot_v1.root'}
def fingerprint(p):
 st=p.stat();h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8388608),b''):h.update(b)
 assert p.stat().st_mtime_ns==st.st_mtime_ns
 return {'path':str(p),'sha256':h.hexdigest(),'bytes':st.st_size,'mtime_ns':st.st_mtime_ns}
provenance={k:fingerprint(p) for k,p in paths.items()};assert provenance['data']['sha256']=='03706f676cf24d3bf0c24c5ab49cf54345369c8cfb80d6f75d3d287e5787e599'
meta=json.loads((base/'artifacts/ppref_x_signed_sweight_all_common.json').read_text());names=[n for n in meta['variables'] if n not in ['Bmass','signal_sWeight']];assert len(names)==67
old={r['variable']:r for r in csv.DictReader((base/'validation/common_scalar_cdf_metrics.csv').open())}
f=R.TFile.Open(str(paths['model']));ws=f.Get('ws_nominal');get=lambda n:ws.var(n).getVal();mean=get('mean1_');sigma=np.sqrt(get('sig1frac1_')*get('sigma11_')**2+(1-get('sig1frac1_'))*get('sigma21_')**2)*get('scale');scale=get('scale');f.Close();windows=[[mean-8*sigma,mean-4*sigma],[mean+4*sigma,mean+8*sigma]]
def read(p,t):
 with uproot.open(p) as f:return f[t].arrays(library='np')
d=read(paths['data'],'ntmix');mc=read(paths['mc'],'ntmix_X3872');sp=read(paths['splot'],'ntmix_X3872_sWeight')
def pre(a):return (a['Bpt']>7.5)&(a['Bpt']<50)&(abs(a['By'])<2.4)&(a['BQvalue']<.15)
sel=pre(d)&(d['Bmass']>3.8)&(d['Bmass']<4.);assert int(sel.sum())==len(sp['Bmass'])
for n in names+['Bmass']:assert np.allclose(d[n][sel],sp[n],atol=1e-7,rtol=0,equal_nan=True),n
sb=pre(d)&(((d['Bmass']>windows[0][0])&(d['Bmass']<windows[0][1]))|((d['Bmass']>windows[1][0])&(d['Bmass']<windows[1][1])))
ms=pre(mc);assert ms.sum()==84381;assert np.unique(mc['pThatreweight']).size==1
samples={'splot':{n:sp[n] for n in names},'mc':{n:mc[n][ms] for n in names},'sideband':{n:d[n][sb] for n in names}};weights={'splot':sp['signal_sWeight'].astype(float),'mc':mc['pThatreweight'][ms].astype(float),'sideband':np.ones(int(sb.sum()))}
config={}
for line in (repo/'plotER/Validation/aux.h').read_text().splitlines():
 q=re.search(r'\{"([^"]+)",\s*"([^"]+)",\s*kNBins,\s*([^,]+),\s*([^,]+),\s*(true|false)',line)
 if q:
  try:config[q[1]]=(float(q[3]),float(q[4]),q[5]=='true',q[2].split(';')[1])
  except ValueError:pass
R.gInterpreter.Declare('#include "'+str(repo/'plotER/Validation/aux.h')+'"')
for q in R.getSignalVars('ntmix_X3872'):config[str(q.expr)]=(float(q.xmin),float(q.xmax),bool(q.absVal),str(q.title).split(';')[1])
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
 rec['cdf_splot_mc']=cdf(*vals['splot'],*vals['mc']);rec['cdf_sideband_mc']=cdf(*vals['sideband'],*vals['mc']);rec['cdf_splot_sideband']=cdf(*vals['splot'],*vals['sideband']);rec['old_cdf_difference']=rec['cdf_splot_mc']-float(old[n]['D_CDF'])
 spH,mcH,sbH=hs['splot'],hs['mc'],hs['sideband']
 for h,color,marker in [(spH,R.kRed+1,24),(sbH,R.kBlue+1,20),(mcH,R.kOrange+7,1)]:h.SetLineColor(color);h.SetMarkerColor(color);h.SetMarkerStyle(marker);h.SetLineWidth(2)
 band=mcH.Clone(n+'_band');band.SetDirectory(0);band.SetFillColorAlpha(R.kOrange+7,.3);band.SetMarkerSize(0)
 c=R.TCanvas('c'+str(j),'',760,650);top=R.TPad('top'+str(j),'',0,.30,1,1);bot=R.TPad('bot'+str(j),'',0,0,1,.30)
 for pad in [top,bot]:pad.SetLeftMargin(.14);pad.SetRightMargin(.04);pad.Draw()
 top.SetBottomMargin(.01);bot.SetTopMargin(.01);bot.SetBottomMargin(.33);top.cd();high=max(h.GetBinContent(i)+h.GetBinError(i) for h in hs.values() for i in range(1,16));low=min(0,min(h.GetBinContent(i)-h.GetBinError(i) for h in hs.values() for i in range(1,16)));spH.SetMaximum(max(high,.1)*1.7);spH.SetMinimum(low*1.12);spH.GetXaxis().SetLabelSize(0);spH.GetYaxis().SetTitle('Normalized entries');spH.Draw('E');band.Draw('E2 SAME');mcH.Draw('HIST SAME');sbH.Draw('E SAME');spH.Draw('E SAME')
 leg=R.TLegend(.57,.76,.94,.91);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.034);leg.AddEntry(mcH,'Prompt X MC','l');leg.AddEntry(spH,'sPlot signal','lep');leg.AddEntry(sbH,'Sideband DATA (raw)','lep');leg.Draw();lab=R.TLatex();lab.SetNDC();lab.SetTextSize(.047);lab.DrawLatex(.18,.87,'ppRef X(3872)');lab.SetTextSize(.031);lab.DrawLatex(.18,.80,'Preselection, no ML cut');lab.DrawLatex(.18,.74,'Sidebands: 4-8#sigma');lab.DrawLatex(.18,.68,'CDF D(sPlot,MC) = %.3f'%rec['cdf_splot_mc'])
 bot.cd();axis=R.TH1D(n+'_ratioaxis','',15,lo,hi);axis.SetDirectory(0);axis.SetMinimum(-2);axis.SetMaximum(5);axis.GetYaxis().SetTitle('DATA / MC');axis.GetYaxis().SetTitleSize(.09);axis.GetYaxis().SetLabelSize(.08);axis.GetYaxis().SetTitleOffset(.7);axis.GetYaxis().SetNdivisions(304);axis.GetXaxis().SetTitle(title);axis.GetXaxis().SetTitleSize(.11);axis.GetXaxis().SetLabelSize(.10);axis.Draw('AXIS');gs=[];clips=0
 for k,color,marker in [('splot',R.kRed+1,24),('sideband',R.kBlue+1,20)]:
  h=hs[k];g=R.TGraphErrors();g.SetName(n+'_'+k+'_ratio');g.SetLineColor(color);g.SetMarkerColor(color);g.SetMarkerStyle(marker)
  for i in range(1,16):
   den=mcH.GetBinContent(i)
   if den<=0:continue
   val=h.GetBinContent(i)/den;err=np.hypot(h.GetBinError(i)/den,h.GetBinContent(i)*mcH.GetBinError(i)/den**2);idx=g.GetN();g.SetPoint(idx,h.GetBinCenter(i),val);g.SetPointError(idx,0,float(err));clips+=int(val< -2 or val>5)
  g.Draw('P SAME');gs.append(g)
 line=R.TLine(lo,1,hi,1);line.SetLineStyle(2);line.Draw();rec['ratio_points_outside_display']=clips
 if clips:lab.SetTextSize(.075);lab.DrawLatex(.18,.85,f'{clips} ratio points outside [-2,5]')
 c.SaveAs(str(out/'figures'/f'{n}.pdf'));root.cd()
 for obj in [*hs.values(),*gs]:obj.Write()
 c.Close();metrics.append(rec)
root.Close();subprocess.run(['pdfunite',*[str(out/'figures'/f'{n}.pdf') for n in names],str(out/'figures/validation_all.pdf')],check=True)
for p in (out/'figures').glob('*.pdf'):assert 'Pages:' in subprocess.check_output(['pdfinfo',str(p)],text=True)
for k,p in paths.items():assert fingerprint(p)==provenance[k]
summary={'status':'PASS','selection':meta['selection'],'sideband_definition':'raw DATA, both sides 4-8 effective sigma, strict open intervals','mean':mean,'effective_sigma':float(sigma),'width_scale':scale,'sideband_windows':windows,'mc_mass_cut':None,'splot_mass_range':[3.8,4.0],'counts':{k:len(w) for k,w in weights.items()},'sumw':{k:float(w.sum()) for k,w in weights.items()},'mc_weight_unique':np.unique(weights['mc']).tolist(),'variables':len(names),'splot_alignment':'all 67 variables and mass matched current DATA entrywise','max_change_from_old_splot_mc_cdf':max(abs(r['old_cdf_difference']) for r in metrics),'input_provenance':provenance,'notes':['MC file hash differs from old snapshot; current constant pThat weighting used.','Existing sPlot and model reused; no refit, no sideband subtraction.','Histogram normalized in display range, nonpositive totals fall back to full support; CDF full finite raw support.','Signed sWeights and negative bins retained; CDF is not a KS p-value.','sPlot scale is near its original upper bound; MC remains single-pThat private sample.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(out/'variable_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n');print(json.dumps(summary,indent=2))
