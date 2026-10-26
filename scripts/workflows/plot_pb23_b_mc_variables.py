#!/usr/bin/env python3
"""PbPb23 B MC variable shapes with the agreed fiducial and MC source weights."""
import argparse, hashlib, json, re, shutil, subprocess
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import uproot
import ROOT as R
from ppref_b_plot_style import labels
parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
out=args.output;assert not out.exists();(out/'metadata').mkdir(parents=True)
repo=Path(__file__).resolve().parents[2]
invpath=repo.parent/'XGBoost/output/internal/code_validation/b_samples_request_20261013/inventory.json'
inv=json.loads(invpath.read_text());template=repo/'plotER/configs/b_inputs/ppref_variables.json';configs=json.loads(template.read_text())
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0);R.gStyle.SetTextFont(42);assert R.gROOT.GetVersion()=='6.32.02'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def stat(p):s=p.stat();return {'size':s.st_size,'mtime_ns':s.st_mtime_ns}
selection='Bpt > 7.5 && abs(By) < 2.4';valid='pthat >= pThats && pThatreweight > 0'
report={'created_UTC':datetime.now(timezone.utc).isoformat(),'status':'running','selection':selection,'MC_validity':valid,'weight':'pThatreweight; retain source5','mass_cut':None,'normalization':'unit visible finite area; sqrt(sumw2) divided by fixed normalization, no normalization covariance','source_SHA_recomputed':False,'inventory_sha256':sha(invpath),'range_template_sha256':sha(template),'script_sha256':sha(Path(__file__)),'channels':{},'visual_review':'pending'}
metrics={};used={};rf=R.TFile(str(out/'metadata/histograms.root'),'RECREATE')
excluded={'RunNo','EvtNo','LumiNo','eventIndex','candidateIndex','Bindex'}
for ch in ['Bu','Bs','Bd']:
 e=next(x for x in inv['inputs'] if x['channel']==ch and x['system']=='PbPb23' and x['kind']=='MC');src=Path(e['path']);before=stat(src);assert before=={k:e[k] for k in before}
 with uproot.open('root://eosuser.cern.ch/'+str(src)) as f:
  t=f[e['tree']];assert t.num_entries==e['entries'];a=t.arrays([n for n in t.keys() if n not in excluded],library='np')
 assert all(np.isfinite(a[n]).all() for n in ['Bpt','By','pthat','pThats','pThatreweight'])
 fid=(a['Bpt'].astype(float)>7.5)&(abs(a['By'].astype(float))<2.4);mask=fid&(a['pthat']>=a['pThats'])&(a['pThatreweight']>0);w=a['pThatreweight'][mask].astype(float);nsel=len(w)
 names=[n for n in configs[ch] if n in a]+sorted(set(a)-set(configs[ch]));used[ch]={};metrics[ch]={};(out/ch/'pdf').mkdir(parents=True)
 df=R.RDataFrame(e['tree'],str(src)).Filter(selection).Filter(valid);count=df.Count();weight=df.Sum('pThatreweight');refs={}
 for n in names:
  x=a[n][mask].astype(float);assert x.ndim==1;good=np.isfinite(x);v=x[good];cfg=dict(configs[ch][n]) if n in configs[ch] else {'bins':15,'range':[0,1],'absolute_display':False,'title':n,'origin':'MC finite support'}
  if n not in configs[ch]:
   if n=='Bmass':cfg.update(range=[5.,6.],bins=100,title={'Bu':'m(J/#psi K^{+}) [GeV/c^{2}]','Bs':'m(J/#psi K^{+}K^{-}) [GeV/c^{2}]','Bd':'m(J/#psi K#pi) [GeV/c^{2}]'}[ch],origin='accepted Bmass 5-6 display')
   elif len(v):
    lo,hi=float(v.min()),float(v.max());pad=(hi-lo)*.025 if hi>lo else .5;cfg['range']=[lo-pad,hi+pad]
  lo,hi=cfg['range'];nb=cfg['bins'];z=abs(v) if cfg['absolute_display'] else v;ww=w[good];edges=np.linspace(lo,hi,nb+1);inside=(z>=lo)&(z<hi)
  indices=np.searchsorted(edges,z[inside],side='right')-1;raw=np.bincount(indices,weights=ww[inside],minlength=nb);sq=np.bincount(indices,weights=ww[inside]**2,minlength=nb)
  under=float(ww[z<lo].sum());over=float(ww[z>=hi].sum());assert np.isclose(raw.sum()+under+over,ww.sum(),rtol=1e-10,atol=1e-7)
  norm=float(raw.sum()) or 1.;rec={'finite':int(good.sum()),'nonfinite':int((~good).sum()),'min':float(v.min()) if len(v) else None,'max':float(v.max()) if len(v) else None,'raw_bins':raw.tolist(),'sumw2_bins':sq.tolist(),'normalization':norm,'underflow_sumw':under,'overflow_sumw':over,'finite_sumw':float(ww.sum()),'inrange_candidates':int(inside.sum()),'inrange_weight_fraction':float(raw.sum()/ww.sum()) if ww.sum()>0 else 0.}
  metrics[ch][n]=rec;used[ch][n]=cfg
  col='display_'+n;dd=df.Define(col,('abs((double)'+n+')') if cfg['absolute_display'] else '(double)'+n).Filter('std::isfinite('+col+')');refs[n]=dd.Histo1D((ch+'_'+n+'_check','',nb,lo,hi),col,'pThatreweight')
 assert int(count.GetValue())==nsel;assert np.isclose(weight.GetValue(),w.sum(),rtol=1e-10)
 for j,n in enumerate(names):
  cfg=used[ch][n];rec=metrics[ch][n];lo,hi=cfg['range'];nb=cfg['bins'];raw=np.array(rec['raw_bins']);sq=np.array(rec['sumw2_bins']);norm=rec['normalization'];hcheck=refs[n].GetValue()
  assert np.allclose(raw,[hcheck.GetBinContent(i) for i in range(1,nb+1)],rtol=1e-9,atol=1e-7),(ch,n,'bins')
  assert np.allclose(sq,[hcheck.GetBinError(i)**2 for i in range(1,nb+1)],rtol=1e-9,atol=1e-6),(ch,n,'errors')
  assert np.allclose([rec['underflow_sumw'],rec['overflow_sumw']],[hcheck.GetBinContent(0),hcheck.GetBinContent(nb+1)],rtol=1e-9,atol=1e-7)
  h=hcheck.Clone(ch+'_'+n+'_raw');h.SetDirectory(0);rf.cd();h.Write();h=h.Clone(ch+'_'+n);h.SetDirectory(0);h.Scale(1/norm);rf.cd();h.Write()
  c=R.TCanvas('c_'+ch+'_'+n,'',760,650);c.SetLeftMargin(.14);c.SetRightMargin(.04);c.SetBottomMargin(.13);c.SetTopMargin(.07)
  h.SetLineColor(R.kOrange+7);h.SetLineWidth(2);h.SetMinimum(0);h.SetMaximum(max(.1,max(h.GetBinContent(i)+h.GetBinError(i) for i in range(1,nb+1)))*1.8);h.GetXaxis().SetTitle(cfg['title']);h.GetYaxis().SetTitle('Normalized entries');h.GetYaxis().SetTitleOffset(1.5);h.Draw('HIST');band=h.Clone('band_'+ch+'_'+n);band.SetDirectory(0);band.SetFillColorAlpha(R.kOrange+7,.3);band.SetMarkerSize(0);band.Draw('E2 SAME');h.Draw('HIST SAME')
  lab=R.TLatex();lab.SetNDC();lab.SetTextFont(42);lab.SetTextSize(.042);lab.DrawLatex(.18,.87,'PbPb 2023 '+labels[ch]+' MC');lab.SetTextSize(.030);lab.DrawLatex(.18,.815,'p_{T} > 7.5 GeV/c, |y| < 2.4');lab.SetTextSize(.026);lab.DrawLatex(.18,.77,f'N = {nsel:,}; in range = {rec["inrange_candidates"]:,}');lab.DrawLatex(.18,.725,f'Weighted outside range: {100*(1-rec["inrange_weight_fraction"]):.2f}%; non-finite: {rec["nonfinite"]}');lab.SetTextSize(.024);lab.DrawLatex(.18,.68,'pThatreweight; valid MC sources; no mass cut')
  if rec['min']==rec['max']:lab.DrawLatex(.18,.64,'Constant branch in selected MC')
  if not raw.sum():lab.DrawLatex(.18,.60,'No candidates in display range')
  c.RedrawAxis();c.Print(str(out/ch/'pdf'/f'{n}.pdf'));c.Close()
 subprocess.run(['pdfunite',*[str(out/ch/'pdf'/f'{n}.pdf') for n in names],str(out/ch/'variables_all.pdf')],check=True)
 assert stat(src)==before
 report['channels'][ch]={'input':e['path'],'tree':e['tree'],**before,'input_candidates':e['entries'],'fiducial_candidates':int(fid.sum()),'selected_candidates':nsel,'sumw':float(w.sum()),'sumw2':float(w@w),'variables':len(names),'ROOT_all_bins_errors_flow':'PASS','source_stat_stable':True}
 print(ch,nsel,len(names),'variables PASS',flush=True)
rf.Close()
with uproot.open(out/'metadata/histograms.root') as f:
 for ch,ms in metrics.items():
  for n,s in ms.items():
   h=f[ch+'_'+n];assert np.allclose(h.values(),np.array(s['raw_bins'])/s['normalization']);assert np.allclose(h.variances(),np.array(s['sumw2_bins'])/s['normalization']**2)
for ch in used:
 for p in list((out/ch/'pdf').glob('*.pdf'))+[out/ch/'variables_all.pdf']:
  info=subprocess.check_output(['pdfinfo',str(p)],text=True);assert int(re.search(r'Pages:\s+(\d+)',info)[1])==(len(used[ch]) if p.name=='variables_all.pdf' else 1)
report['status']='numerical_PASS';report['pdf_page_counts']='PASS';report['persisted_bins_errors']='PASS'
for n,x in [('provenance.json',report),('variable_metrics.json',metrics),('variables.json',used)]:dump(out/'metadata'/n,x)
for p in [Path(__file__),invpath]:shutil.copy2(p,out/'metadata'/p.name)
(out/'README.md').write_text('# Bvar_pb23_MC\n\nPbPb23 Bu/Bs/Bd MC，基础fiducial：Bpt > 7.5 GeV、abs(By) < 2.4。无2D、mass、中心度或额外channel pre-cut。沿用MC来源有效性pthat >= pThats、pThatreweight > 0，以pThatreweight加权，保留来源5。\n\n'+''.join(f'- [{ch}合订本]({ch}/variables_all.pdf) · [单变量PDF]({ch}/pdf/)\n' for ch in used)+'\n沿用已有Variables_ppRef变量范围、分箱及MC橙色曲线/误差带，新增MC字段按完整有限支持区间显示，Bmass显示5–6 GeV、10 MeV/bin。排除运行/事件/候选索引，其余标量含MC诊断字段全部绘制。每图按可见有限范围单位面积归一化，误差sqrt(sumw2)/固定归一化常数，不包含归一化协方差。N和in range为未加权候选数；范围外比例按权重计算。非有限值仅从对应变量排除；恒定占位字段不推定物理有效性。\n\n[输入与核验](metadata/provenance.json) · [变量统计](metadata/variable_metrics.json) · [范围](metadata/variables.json)。ROOT独立核对所有变量bins、sumw2和flow；源size/mtime读前后及inventory一致，未重算大文件SHA。仅样式重绘可复用histograms.root；输入、selection、权重、变量或分箱改变必须重算。\n')
print('ALL CHANNELS NUMERICALLY VALIDATED',flush=True)
