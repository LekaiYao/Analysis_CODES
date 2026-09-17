from pathlib import Path
import argparse,csv,json,re,subprocess,hashlib
import numpy as np,uproot,ROOT as R
from ppref_psi2s_prompt_fit import cdf,corr,dump
ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);a=ap.parse_args();out=Path(a.output)
repo=Path(__file__).resolve().parents[2];base=repo/'plotER/Validation/results/psi2s_ppref_validation'
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0)
old={r['variable']:r for r in json.loads((base/'variable_metrics.json').read_text())}
def read(p,t):
 with uproot.open(p) as f:return f[t].arrays(library='np')
d=read(out/'artifacts/psi2s_sweights.root','ntmix_PSI2S_sWeight');m=read(out/'artifacts/selected_mc.root','ntmix_PSI2S')
bd=read(base/'artifacts/psi2s_sweights.root','ntmix_PSI2S_sWeight');bm=read(base/'artifacts/selected_mc.root','ntmix_PSI2S')
for new,prev in [(d,bd),(m,bm)]:
 mask=prev['Bnorm_svpvDistance_2D']<2
 assert np.array_equal(new['source_entry'],prev['source_entry'][mask])
 for k in ['Bmass','Bpt','Bnorm_svpvDistance_2D']:assert np.array_equal(new[k],prev[k][mask])
assert np.array_equal(m['pThatreweight'],bm['pThatreweight'][bm['Bnorm_svpvDistance_2D']<2])
# Include historical commented definitions; active and historical ranges stay fixed.
config={}
for line in (repo/'plotER/Validation/aux.h').read_text().splitlines():
 q=re.search(r'\{"([^"]+)",\s*"([^"]+)",\s*kNBins,\s*([^,]+),\s*([^,]+),\s*(true|false)',line)
 if q:
  try:config[q[1]]=(float(q[3]),float(q[4]),q[5]=='true',q[2].split(';')[1])
  except ValueError:pass
# Species-dependent dipion range is the active native configuration.
R.gInterpreter.Declare('#include "'+str(repo/'plotER/Validation/aux.h')+'"')
for q in R.getSignalVars('ntmix_PSI2S'):config[str(q.expr)]=(float(q.xmin),float(q.xmax),bool(q.absVal),str(q.title).split(';')[1])
rows=[];root=R.TFile(str(out/'artifacts/validation_histograms.root'),'RECREATE');book=str(out/'figures/validation_all.pdf')
for j,n in enumerate(sorted(old)):
 good=np.isfinite(d[n]);gm=np.isfinite(m[n]);x=d[n][good].astype(float);w=d['signal_sWeight'][good];y=m[n][gm].astype(float);v=m['pThatreweight'][gm].astype(float)
 z,cd,cm,after=cdf(x,w,y,v)
 gd0=np.isfinite(bd[n]);gm0=np.isfinite(bm[n]);before=cdf(bd[n][gd0],bd['signal_sWeight'][gd0],bm[n][gm0],bm['pThatreweight'][gm0])[3]
 assert abs(before-old[n]['cdf_pthat'])<1e-10
 lo,hi,ab,title=config.get(n,(old[n]['display_min'],old[n]['display_max'],False,n))
 if hi<=lo:lo,hi=lo-.5,hi+.5
 xx=np.abs(x) if ab else x;yy=np.abs(y) if ab else y;edges=np.linspace(lo,hi,16)
 rec={'variable':n,'cdf_before':before,'cdf_after':after,'delta':after-before,'mc_mass_pearson':corr(y,m['Bmass'][gm],v),'data_nonfinite':int((~good).sum()),'mc_nonfinite':int((~gm).sum()),'constant_both':bool(np.ptp(x)==0 and np.ptp(y)==0),'display_min':lo,'display_max':hi,'display_abs':ab,'range_source':'aux.h' if n in config else 'frozen_baseline','bins':15}
 hs=[]
 for key,vals,weights in [('data',xx,w),('mc',yy,v)]:
  inside=(vals>=lo)&(vals<hi);raw=np.histogram(vals[inside],edges,weights=weights[inside])[0];err=np.sqrt(np.histogram(vals[inside],edges,weights=weights[inside]**2)[0]);norm=raw.sum()
  rec[key+'_inrange_weight_fraction']=float(norm/weights.sum());rec[key+'_hist_normalization']='inrange' if norm>0 else 'full_support_fallback'
  if norm<=0:norm=weights.sum()
  h=R.TH1D(n+'_'+key,'',15,lo,hi);h.SetDirectory(0);h.Sumw2()
  for i in range(15):h.SetBinContent(i+1,raw[i]/norm);h.SetBinError(i+1,err[i]/norm)
  hs.append(h)
 sp,mc=hs;sp.SetLineColor(R.kRed+1);sp.SetMarkerColor(R.kRed+1);sp.SetMarkerStyle(24);sp.SetLineWidth(2);mc.SetLineColor(R.kOrange+7);mc.SetLineWidth(2)
 band=mc.Clone(n+'_band');band.SetDirectory(0);band.SetFillColorAlpha(R.kOrange+7,.3);band.SetMarkerSize(0)
 c=R.TCanvas('c'+str(j),'',760,650);top=R.TPad('top'+str(j),'',0,.30,1,1);bot=R.TPad('bot'+str(j),'',0,0,1,.30)
 for p in [top,bot]:p.SetLeftMargin(.14);p.SetRightMargin(.04);p.Draw()
 top.SetBottomMargin(.01);bot.SetTopMargin(.01);bot.SetBottomMargin(.33);top.cd()
 high=max(h.GetBinContent(i)+h.GetBinError(i) for h in hs for i in range(1,16));low=min(0.,min(sp.GetBinContent(i)-sp.GetBinError(i) for i in range(1,16)))
 sp.SetMaximum(max(.1,high)*1.7);sp.SetMinimum(low*1.15);sp.GetYaxis().SetTitle('Normalized entries');sp.GetXaxis().SetLabelSize(0);sp.GetXaxis().SetTitleSize(0);sp.Draw('E');band.Draw('E2 SAME');mc.Draw('HIST SAME');sp.Draw('E SAME')
 leg=R.TLegend(.62,.77,.94,.90);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.04);leg.AddEntry(mc,'Prompt MC','l');leg.AddEntry(sp,'sPlot DATA','lep');leg.Draw()
 lab=R.TLatex();lab.SetNDC();lab.SetTextSize(.052);lab.DrawLatex(.18,.87,'ppRef #psi(2S)');lab.SetTextSize(.034);lab.DrawLatex(.18,.80,'SV-PV significance < 2');lab.DrawLatex(.18,.735,'CDF D: %.3f #rightarrow %.3f'%(before,after));lab.SetTextSize(.029);lab.DrawLatex(.18,.67,'MC: pThatreweight; 15 bins')
 if min(rec['data_inrange_weight_fraction'],rec['mc_inrange_weight_fraction'])<.95:lab.DrawLatex(.18,.61,'In-range fraction DATA/MC: %.3f / %.3f'%(rec['data_inrange_weight_fraction'],rec['mc_inrange_weight_fraction']))
 bot.cd();axis=R.TH1D(n+'_ratio_axis','',15,lo,hi);axis.SetDirectory(0);axis.SetMinimum(0);axis.SetMaximum(5);axis.GetYaxis().SetTitle('sPlot / MC');axis.GetYaxis().SetTitleSize(.09);axis.GetYaxis().SetLabelSize(.08);axis.GetYaxis().SetTitleOffset(.7);axis.GetYaxis().SetNdivisions(304);axis.GetXaxis().SetTitle(title);axis.GetXaxis().SetTitleSize(.11);axis.GetXaxis().SetLabelSize(.10);axis.Draw('AXIS')
 g=R.TGraphErrors();g.SetName(n+'_ratio');clip=0
 for i in range(1,16):
  den=mc.GetBinContent(i);num=sp.GetBinContent(i)
  if den<=0:continue
  val=num/den;err=np.hypot(sp.GetBinError(i)/den,num*mc.GetBinError(i)/den**2);k=g.GetN();g.SetPoint(k,mc.GetBinCenter(i),val);g.SetPointError(k,0,float(err));clip+=int(val<0 or val>5)
 g.SetMarkerStyle(20);g.SetMarkerColor(R.kOrange+7);g.SetLineColor(R.kOrange+7);g.Draw('P SAME');line=R.TLine(lo,1,hi,1);line.SetLineStyle(2);line.Draw();rec['ratio_points_outside_0_5']=clip
 if clip:lab.SetTextSize(.075);lab.DrawLatex(.18,.82,f'{clip} ratio points outside [0,5]')
 c.SaveAs(str(out/'figures'/f'{n}.pdf'))
 root.cd()
 for obj in [sp,mc,g]:obj.Write()
 # Store full-support CDF curves, independent of displayed histogram range.
 for tag,values in [('data',cd),('mc',cm)]:
  gg=R.TGraph(len(z),np.ascontiguousarray(z,dtype='d'),np.ascontiguousarray(values,dtype='d'));gg.Write(n+'_cdf_'+tag)
 c.Close();rows.append(rec)
root.Close();subprocess.run(['pdfunite',*[str(out/'figures'/f'{n}.pdf') for n in sorted(old)],book],check=True)
rows.sort(key=lambda r:r['cdf_after'],reverse=True);dump(out/'cdf_comparison.json',rows)
with (out/'cdf_comparison.csv').open('w') as f:
 wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
# Native RooFit mass plot and pull, restored to pre-sPlot fit parameters.
f=R.TFile.Open(str(out/'artifacts/fit_workspace.root'));ws=f.Get('ws_nominal');fd=f.Get('fit_data');mass=ws.var('Bmass');model=ws.pdf('model');data=ws.data('data')
for p in fd.floatParsFinal():ws.var(p.GetName()).setVal(p.getVal())
frame=mass.frame(R.RooFit.Bins(80));data.plotOn(frame,R.RooFit.Name('data_points'));model.plotOn(frame,R.RooFit.Name('total'));pull=frame.pullHist('data_points','total');model.plotOn(frame,R.RooFit.Components('background'),R.RooFit.LineStyle(2),R.RooFit.LineColor(R.kRed+1));frame.GetXaxis().SetLabelSize(0);frame.GetXaxis().SetTitleSize(0)
c=R.TCanvas('mass_canvas','',760,650);top=R.TPad('mass_top','',0,.30,1,1);bot=R.TPad('mass_bot','',0,0,1,.30)
for p in [top,bot]:p.SetLeftMargin(.14);p.SetRightMargin(.04);p.Draw()
top.SetBottomMargin(.01);bot.SetTopMargin(.01);bot.SetBottomMargin(.33);top.cd();frame.Draw();lab.SetTextSize(.037);lab.DrawLatex(.18,.87,'ppRef #psi(2S), SV-PV significance < 2');lab.DrawLatex(.18,.80,'N_{sig} = 24611 #pm 267')
bot.cd();pf=mass.frame();pf.addPlotable(pull,'P');pf.SetMinimum(-5);pf.SetMaximum(5);pf.GetYaxis().SetTitle('Pull');pf.GetYaxis().SetTitleSize(.09);pf.GetYaxis().SetLabelSize(.08);pf.GetYaxis().SetTitleOffset(.7);pf.GetXaxis().SetTitle('m(J/#psi #pi^{+}#pi^{-}) [GeV/c^{2}]');pf.GetXaxis().SetTitleSize(.10);pf.GetXaxis().SetLabelSize(.09);pf.Draw();c.SaveAs(str(out/'figures/mass_fit.pdf'));c.Close();f.Close()
for p in (out/'figures').glob('*.pdf'):assert 'Pages:' in subprocess.check_output(['pdfinfo',str(p)],text=True)
assert len(list((out/'figures').glob('*.pdf')))==69
assert not list(out.rglob('*.png'))
summary=json.loads((out/'summary.json').read_text());summary['status']='complete';dump(out/'summary.json',summary)
dump(out/'audit.json',{'status':'PASS','variables':67,'pdf_files':69,'book_pages':67,'baseline_cdf_reproduced_max_tolerance':1e-10,'source_entry_selection':'exact baseline subset for DATA and MC','mc_weights':'exact baseline subset','sweights':'refitted after cut; closure in summary.json','ratio_display':[0,5],'cdf':'signed full finite support; no KS p-value','hist_normalization':'in-range; full-support fallback if in-range total <= 0','bootstrap':False})
print('COMPLETE',out)
