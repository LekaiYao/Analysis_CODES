# Historical study: frozen tags/FOM contract; not the current pb23 workflow.
from pathlib import Path
import json,csv,re,hashlib,subprocess,argparse
import numpy as np,uproot,ROOT as R
ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--cache-source',type=Path);ap.add_argument('--no-dr-cut',action='store_true');ap.add_argument('--variables',nargs='+');ap.add_argument('--trk2-dr-only',action='store_true');args=ap.parse_args();assert not (args.no_dr_cut and args.trk2_dr_only);repo=Path.cwd();out=Path(args.output);out.mkdir(parents=True,exist_ok=False);(out/'pdf').mkdir();(out/'artifacts').mkdir()
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0)
src=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/PbPb23');paths={'data':src/'flat_ntmix_PbPb23_DATA.root','mc':src/'flat_ntmix_PbPb23_MC_X3872.root'}
def fingerprint(p):
 st=p.stat();return {'path':str(p),'bytes':st.st_size,'mtime_ns':st.st_mtime_ns}
provenance={k:fingerprint(p) for k,p in paths.items()}
requested='Bmu1pt Bmu2pt Btrk1Pt Btrk2Pt Btrk1dR Btrk2dR BtrkPtimb Btktkpt Bchi2Prob BtktkvProb Balpha Bcos_dtheta'.split()
allowed=set(requested+'Bmu1eta Bmu2eta Btrk1Eta Btrk2Eta Bujeta Btktketa BujvProb Bnorm_svpvDistance_2D BsvpvDistance_2D Bdtheta Bujpt Btrk1Dxy1 Btrk2Dxy1 Bnorm_trk1Dxy Bnorm_trk2Dxy'.split())
if args.variables:requested=args.variables
missing=[];names=list(requested);assert names and len(names)==len(set(names)) and set(names)<=allowed
windows=[[3.80,3.85],[3.89,3.94]]
selection='Bpt > 15 && Bpt < 50 && abs(By) < 1.6 && BQvalue < 0.15'+('' if args.no_dr_cut else (' && Btrk2dR < 0.35' if args.trk2_dr_only else ' && Btrk1dR < 0.35 && Btrk2dR < 0.35'))
cache_names=list(dict.fromkeys(names+['Bmass','Bpt','By','BQvalue']))
template_path=repo/'plotER/validation_archive/archive/previous_validation/results/pb23_x_mc_sideband/variable_metrics.json'
template={r['variable']:r for r in json.loads(template_path.read_text())};assert set(names)<=set(template)
samples={};weights={};cutflow={}
with uproot.open(paths['mc']) as f:
 source_weights=f['ntmix_X3872']['pThatreweight'].array(library='np');assert np.isfinite(source_weights).all() and (source_weights>0).all()
 excluded_weight=float(source_weights.max());source_weight_values=np.unique(source_weights).tolist()
if args.cache_source:
 prior=json.loads((args.cache_source/'summary.json').read_text())
 assert prior['selection']==selection
 assert prior['sideband_windows']==windows and set(names)<=set(prior['plotted_variables']) and prior['mc_mass_cut'] is None
 assert prior['input_provenance']==provenance
 for k in ['mc','sideband']:
  cache=args.cache_source/'artifacts'/f'{k}_selected.npz'
  hashes=json.loads((args.cache_source/'delivery_hashes.json').read_text());assert hashlib.sha256(cache.read_bytes()).hexdigest()==hashes[str(cache.relative_to(args.cache_source))]
  arrays=dict(np.load(cache));weights[k]=arrays.pop('event_weight');samples[k]={n:arrays[n] for n in names};cutflow[k]=prior['cutflow'][k].copy()
else:
 for key,tree in [('mc','ntmix_X3872'),('data','ntmix')]:
  chunks={n:[] for n in cache_names};wc=[];total=0;precount=0;passed=0
  with uproot.open(paths[key]) as f:
   tr=f[tree];assert set(names).issubset(tr.keys());read=sorted(set(names+['Bmass','Bpt','By','BQvalue']+(['pThatreweight'] if key=='mc' else [])))
   for a in tr.iterate(read,step_size=100000,library='np'):
    pre=(a['Bpt']>15)&(a['Bpt']<50)&(abs(a['By'])<1.6)&(a['BQvalue']<.15)
    if not args.no_dr_cut:
     pre &= a['Btrk2dR']<.35
     if not args.trk2_dr_only:pre &= a['Btrk1dR']<.35
    total+=len(pre);precount+=int(pre.sum());mask=pre.copy()
    if key=='data':mask &= ((a['Bmass']>3.80)&(a['Bmass']<3.85))|((a['Bmass']>3.89)&(a['Bmass']<3.94))
    passed+=int(mask.sum())
    for n in cache_names:chunks[n].append(a[n][mask])
    w=a['pThatreweight'][mask].astype(float) if key=='mc' else np.ones(int(mask.sum()));assert np.isfinite(w).all() and (w>0).all();wc.append(w)
  label='sideband' if key=='data' else 'mc';samples[label]={n:np.concatenate(v) for n,v in chunks.items()};weights[label]=np.concatenate(wc);cutflow[label]={'total':total,'preselection':precount,'selected':passed};print('SELECTED',label,cutflow[label],flush=True)
mc_keep=weights['mc']!=excluded_weight
if args.cache_source and prior.get('mc_excluded_weight')==excluded_weight:
 assert mc_keep.all()
 excluded_entries=prior['mc_excluded_entries']
 assert cutflow['mc']['before_pthat5_veto']-cutflow['mc']['selected']==excluded_entries
else:
 excluded_entries=int((~mc_keep).sum());assert excluded_entries>0
 samples['mc']={n:x[mc_keep] for n,x in samples['mc'].items()};weights['mc']=weights['mc'][mc_keep]
 cutflow['mc']['before_pthat5_veto']=cutflow['mc']['selected'];cutflow['mc']['pthat5_excluded']=excluded_entries;cutflow['mc']['selected']=int(mc_keep.sum())
assert np.all(weights['mc']<excluded_weight)
# Reusable selected arrays: invalidated by source stat, exact selection, weights or variable list changes.
for k in samples:np.savez_compressed(out/'artifacts'/f'{k}_selected.npz',**samples[k],event_weight=weights[k])
old={}
for n in names:
 vals=np.r_[samples['mc'][n],samples['sideband'][n]];vals=vals[np.isfinite(vals)];old[n]={'support_min':float(vals.min()),'support_max':float(vals.max())}
config={}
for line in (repo/'plotER/validation_archive/aux.h').read_text().splitlines():
 q=re.search(r'\{"([^"]+)",\s*"([^"]+)",\s*kNBins,\s*([^,]+),\s*([^,]+),\s*(true|false)',line)
 if q:
  try:config[q[1]]=(float(q[3]),float(q[4]),q[5]=='true',q[2].split(';')[1])
  except ValueError:pass
R.gInterpreter.Declare('#include "'+str(repo/'plotER/validation_archive/aux.h')+'"')
for q in R.getSignalVars('ntmix_X3872'):config[str(q.expr)]=(float(q.xmin),float(q.xmax),bool(q.absVal),str(q.title).split(';')[1])
if 'Btktkmass' in names:
 lo=np.floor(old['Btktkmass']['support_min']/.05)*.05;hi=(np.floor(old['Btktkmass']['support_max']/.05)+1)*.05
 config['Btktkmass']=(float(lo),float(hi),False,'m_{#pi#pi} [GeV/c^{2}]')
def cdf(x,w,y,v):
 z=np.union1d(x,y);ix=np.argsort(x,kind='stable');iy=np.argsort(y,kind='stable');a=np.r_[0,np.cumsum(w[ix])][np.searchsorted(x[ix],z,side='right')]/w.sum();b=np.r_[0,np.cumsum(v[iy])][np.searchsorted(y[iy],z,side='right')]/v.sum();return float(np.max(abs(a-b)))
metrics=[];root=R.TFile(str(out/'artifacts/histograms.root'),'RECREATE')
for j,n in enumerate(names):
 lo,hi=template[n]['range'];ab=template[n]['absolute_display'];title=config.get(n,(lo,hi,ab,n))[3]
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
 leg=R.TLegend(.57,.76,.94,.91);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.034);leg.AddEntry(mcH,'Prompt X MC (no pThat5)','l');leg.AddEntry(sbH,'Sideband DATA (raw)','lep');leg.Draw();lab=R.TLatex();lab.SetNDC();lab.SetTextSize(.047);lab.DrawLatex(.18,.87,'PbPb23 X(3872)');lab.SetTextSize(.031);lab.DrawLatex(.18,.80,('No dR or ML cut' if args.no_dr_cut else ('trk2 dR < 0.35, no ML cut' if args.trk2_dr_only else 'Preselection, no ML cut')));lab.DrawLatex(.18,.74,'SB: (3.80,3.85) + (3.89,3.94)');lab.DrawLatex(.18,.68,'CDF D(SB,MC) = %.3f'%rec['cdf_sideband_mc'])
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
 c.SaveAs(str(out/'pdf'/f'{n}.pdf'));root.cd()
 for obj in [*hs.values(),*gs]:obj.Write()
 c.Close();metrics.append(rec)
root.Close();subprocess.run(['pdfunite',*[str(out/'pdf'/f'{n}.pdf') for n in names],str(out/'validation_all.pdf')],check=True)
for p in [*(out/'pdf').glob('*.pdf'),out/'validation_all.pdf']:assert 'Pages:' in subprocess.check_output(['pdfinfo',str(p)],text=True)
for k,p in paths.items():assert fingerprint(p)==provenance[k]
summary={'status':'PASS','selection':selection,'dr_cut_applied':not args.no_dr_cut,'trk2_dr_only':args.trk2_dr_only,'sideband_definition':'raw DATA, strict open mass intervals','sideband_windows':windows,'mc_mass_cut':None,'counts':{k:len(w) for k,w in weights.items()},'cutflow':cutflow,'sumw':{k:float(w.sum()) for k,w in weights.items()},'mc_weight_unique':np.unique(weights['mc']).tolist(),'mc_excluded_weight':excluded_weight,'mc_excluded_entries':excluded_entries,'source_mc_weight_values':source_weight_values,'variables':len(names),'requested_variables':requested,'plotted_variables':names,'missing_variables':missing,'variable_scope':'Explicit variable list recorded in requested_variables','template_metrics':str(template_path),'input_provenance':provenance,'notes':['Prompt X MC pThatreweight, excluding maximum-weight pthat5 category; raw sideband DATA unit weight.',('No dR, ML or displacement cut; no fitting or subtraction.' if args.no_dr_cut else 'No ML or displacement cut; dR upper bounds only, no additional lower bounds; no fitting or subtraction.'),'Histogram normalized in visible range; full-support fallback for empty visible range recorded per variable.','CDF full finite raw support; separation diagnostic, not MC-to-signal validation.','Input file size/mtime checked before/after; no claim of full source SHA256.']}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(out/'variable_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n');print(json.dumps(summary,indent=2))

import shutil
for p in [Path(__file__),repo/'plotER/validation_archive/aux.h',template_path]:shutil.copy2(p,out/'artifacts'/p.name)
(out/'README.md').write_text('# PbPb23 X prompt MC vs raw sideband\n\n[合订本（'+str(len(names))+'页）](validation_all.pdf) · [单页PDF](pdf/) · [筛选与计数](summary.json)\n\n本次指定的'+str(len(names))+'项变量全部完成。\n\n共同筛选：'+selection+'。无ML及位移筛选；是否有dR cut以此明确表达式为准。DATA侧带为(3.80,3.85)或(3.89,3.94)，严格开区间；MC不加质量窗，按pThatreweight填图，但排除其最大权重一档（pthat5）；DATA单位权重。\n\n沿用已接受的15 bins、逐变量范围、可见范围归一化、MC误差带及Sideband/MC比值[0,5]。范围外权重比例逐变量记录；CDF使用完整有限原始取值。MC与背景侧带的形状差异是诊断对照，不能据此直接判定MC对真实signal的描述质量。\n\nartifacts内保留筛选后数组与冻结绘图源码；缓存只适用于相同源文件stat、筛选、权重和变量定义。\n')
(out/'delivery_hashes.json').write_text(json.dumps({str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob('*') if p.is_file()},indent=2)+'\n')
