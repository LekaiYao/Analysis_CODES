#!/usr/bin/env python3
"""Bu ppRef raw DATA versus current MC; other channels await current ppRef MC."""
from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,shutil,subprocess,re
import numpy as np
import uproot
import ROOT as R
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0);R.gStyle.SetOptTitle(0);R.gStyle.SetTextFont(42)
assert R.gROOT.GetVersion()=='6.32.02'
repo=Path(__file__).resolve().parents[2];out=Path('/tmp/leyao/Bvar_ppRef_dataMC_2D4');assert not out.exists();(out/'Bu/pdf').mkdir(parents=True);(out/'metadata').mkdir()
previous=repo/'plotER/configs/b_inputs';variable_order=json.loads((previous/'variable_order.json').read_text());configs=json.loads((previous/'ppref_variables.json').read_text())['Bu'];names=variable_order['Bu']
inventory=json.loads((repo.parent/'XGBoost/output/internal/code_validation/b_samples_request_20261013/inventory.json').read_text())
source={key:next(x for x in inventory['inputs'] if x['system']=='ppRef' and x['channel']=='Bu' and x['kind']==key) for key in ['DATA','MC']}
selection='Bpt > 7.5 && abs(By) < 2.4 && Bnorm_svpvDistance_2D > 4';mc_valid='pthat >= pThats && pThatreweight > 0'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
report={'created_UTC':datetime.now(timezone.utc).isoformat(),'status':'partial_inputs','channels':{'Bu':'running','Bs':'awaiting current ppRef MC','Bd':'awaiting current ppRef MC'},'fiducial':'Bpt > 7.5 && abs(By) < 2.4','selection':selection,'mass_cut':None,'MC_validity':mc_valid,'MC_weight':'pThatreweight, retain source5','DATA':'raw inclusive candidates; signal + background; no subtraction','normalization':'visible finite range independently for DATA and MC; empty visible distribution uses full finite sumw only to define empty histogram','error':'sqrt(sumw2)/fixed normalization; ratio independent DATA and MC propagation, no normalization covariance','template':str(repo/'scripts/workflows/ppref_psi2s_prompt_plots.py'),'template_sha256':sha(repo/'scripts/workflows/ppref_psi2s_prompt_plots.py'),'range_template':str(previous/'ppref_variables.json'),'range_template_sha256':sha(previous/'ppref_variables.json'),'script_sha256':sha(Path(__file__)),'samples':{}}
histos={};summaries={}
for key,record in source.items():
 src=Path(record['path']);before=src.stat();assert (before.st_size,before.st_mtime_ns)==(record['size'],record['mtime_ns'])
 hists={n:np.zeros(15) for n in names};sumw2={n:np.zeros(15) for n in names};stats={n:dict(finite=0,nonfinite=0,sumw=0.,underflow_sumw=0.,overflow_sumw=0.,min=None,max=None) for n in names}
 ninput=nfid=nsel=0;sw=sw2=0.
 with uproot.open(src) as f:
  t=f[record['tree']];assert set(names)<=set(t.keys())
  columns=list(dict.fromkeys(names+(['pthat','pThats','pThatreweight'] if key=='MC' else [])))
  for a in t.iterate(columns,step_size='80 MB',library='np'):
   pt=a['Bpt'].astype(float);y=a['By'].astype(float);d=a['Bnorm_svpvDistance_2D'].astype(float)
   assert np.isfinite(pt).all() and np.isfinite(y).all() and np.isfinite(d).all()
   mask=(pt>7.5)&(abs(y)<2.4)&(d>4);ninput+=len(pt);nfid+=int(mask.sum())
   if key=='MC':
    assert all(np.isfinite(a[n]).all() for n in ['pthat','pThats','pThatreweight']);mask&=(a['pthat']>=a['pThats'])&(a['pThatreweight']>0)
   w=a['pThatreweight'][mask].astype(float) if key=='MC' else np.ones(int(mask.sum()));nsel+=len(w);sw+=float(w.sum());sw2+=float(w@w)
   for n in names:
    x=a[n][mask].astype(float);good=np.isfinite(x);s=stats[n];s['nonfinite']+=int((~good).sum());s['finite']+=int(good.sum());x=x[good];v=w[good];s['sumw']+=float(v.sum())
    if len(x):s['min']=float(x.min()) if s['min'] is None else min(s['min'],float(x.min()));s['max']=float(x.max()) if s['max'] is None else max(s['max'],float(x.max()))
    cfg=configs[n];lo,hi=cfg['range'];x=abs(x) if cfg['absolute_display'] else x;inside=(x>=lo)&(x<hi);edges=np.linspace(lo,hi,16)
    indices=np.searchsorted(edges,x[inside],side='right')-1;hists[n]+=np.bincount(indices,weights=v[inside],minlength=15);sumw2[n]+=np.bincount(indices,weights=v[inside]**2,minlength=15);s['underflow_sumw']+=float(v[x<lo].sum());s['overflow_sumw']+=float(v[x>=hi].sum())
 assert ninput==record['entries']
 if key=='DATA':assert nsel==273455
 # Independent ROOT selection and weighted total, plus three principal distributions.
 df=R.RDataFrame(record['tree'],str(src)).Filter(selection)
 if key=='MC':df=df.Filter(mc_valid)
 count=df.Count();weight=df.Sum('pThatreweight') if key=='MC' else None;refs={}
 for n in ['Bpt','By','Bnorm_svpvDistance_2D']:
  cfg=configs[n];display='abs_'+n if cfg['absolute_display'] else n;dd=df.Define(display,'abs('+n+')') if display!=n else df
  refs[n]=dd.Histo1D((key+'_'+n+'_check','',15,*cfg['range']),display,'pThatreweight') if key=='MC' else dd.Histo1D((key+'_'+n+'_check','',15,*cfg['range']),display)
 assert int(count.GetValue())==nsel
 if weight:assert np.isclose(weight.GetValue(),sw,rtol=1e-10)
 for n,hp in refs.items():
  h=hp.GetValue();assert np.allclose(hists[n],[h.GetBinContent(i) for i in range(1,16)],rtol=1e-9,atol=1e-7);assert np.allclose(sumw2[n],[h.GetBinError(i)**2 for i in range(1,16)],rtol=1e-9,atol=1e-6)
 for n,s in stats.items():assert s['finite']+s['nonfinite']==nsel;assert np.isclose(sum(hists[n])+s['underflow_sumw']+s['overflow_sumw'],s['sumw'],rtol=1e-10,atol=1e-7)
 after=src.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
 report['samples'][key]={'path':str(src),'tree':record['tree'],'size':before.st_size,'mtime_ns':before.st_mtime_ns,'input_candidates':ninput,'fiducial_and_2D4_candidates':nfid,'selected_candidates':nsel,'sumw':sw,'sumw2':sw2,'ROOT_selection_and_principal_bins':'PASS','source_stat_stable':True}
 histos[key]=(hists,sumw2);summaries[key]=stats
 print(key,nsel,'sumw',sw,'PASS',flush=True)
root=R.TFile(str(out/'metadata/histograms.root'),'RECREATE');metrics={}
for j,n in enumerate(names):
 cfg=configs[n];lo,hi=cfg['range'];hs={};rec={'config':cfg,'samples':{},'ratio_points':[]}
 for key in ['DATA','MC']:
  raw,squares=histos[key][0][n],histos[key][1][n];s=summaries[key][n];norm=float(raw.sum());fallback=norm<=0
  if fallback:norm=s['sumw'] if s['sumw']>0 else 1.
  h=R.TH1D(n+'_'+key,'',15,lo,hi);h.SetDirectory(0);h.Sumw2();r=h.Clone(n+'_'+key+'_raw');r.SetDirectory(0)
  for i in range(15):h.SetBinContent(i+1,raw[i]/norm);h.SetBinError(i+1,np.sqrt(squares[i])/norm);r.SetBinContent(i+1,raw[i]);r.SetBinError(i+1,np.sqrt(squares[i]))
  root.cd();r.Write();h.Write();hs[key]=h;rec['samples'][key]={**s,'inrange_fraction':float(raw.sum()/s['sumw']) if s['sumw']>0 else None,'normalization':norm,'empty_visible_range':fallback,'raw_bins':raw.tolist(),'sumw2_bins':squares.tolist()}
 data,mc=hs['DATA'],hs['MC'];data.SetLineColor(R.kBlue+1);data.SetMarkerColor(R.kBlue+1);data.SetMarkerStyle(20);data.SetMarkerSize(.65);data.SetLineWidth(2);mc.SetLineColor(R.kOrange+7);mc.SetLineWidth(2)
 band=mc.Clone(n+'_band');band.SetDirectory(0);band.SetFillColorAlpha(R.kOrange+7,.3);band.SetMarkerSize(0)
 c=R.TCanvas('c'+str(j),'',760,650);top=R.TPad('top'+str(j),'',0,.30,1,1);bot=R.TPad('bot'+str(j),'',0,0,1,.30)
 for pad in [top,bot]:pad.SetLeftMargin(.14);pad.SetRightMargin(.04);pad.Draw()
 top.SetBottomMargin(.01);bot.SetTopMargin(.01);bot.SetBottomMargin(.33);top.cd()
 high=max(h.GetBinContent(i)+h.GetBinError(i) for h in hs.values() for i in range(1,16));data.SetMaximum(max(high,.1)*1.9);data.SetMinimum(0);data.GetXaxis().SetLabelSize(0);data.GetYaxis().SetTitle('Normalized entries');data.Draw('E');band.Draw('E2 SAME');mc.Draw('HIST SAME');data.Draw('E SAME')
 leg=R.TLegend(.66,.78,.95,.91);leg.SetBorderSize(0);leg.SetFillStyle(0);leg.SetTextSize(.035);leg.AddEntry(data,'Raw DATA','lep');leg.AddEntry(mc,'MC (weighted)','l');leg.Draw()
 lab=R.TLatex();lab.SetNDC();lab.SetTextFont(42);lab.SetTextSize(.043);lab.DrawLatex(.18,.87,'ppRef B^{+}');lab.SetTextSize(.029)
 lab.DrawLatex(.18,.81,'p_{T} > 7.5 GeV/c, |y| < 2.4');lab.DrawLatex(.18,.76,'Bnorm_svpvDistance_2D > 4');lab.SetTextSize(.026);lab.DrawLatex(.18,.705,'No mass cut; raw DATA includes background')
 fd=rec['samples']['DATA']['inrange_fraction'];fm=rec['samples']['MC']['inrange_fraction'];lab.DrawLatex(.18,.655,f'In-range DATA/MC: {100*(fd or 0):.1f}% / {100*(fm or 0):.1f}%')
 extras=[]
 for key in ['DATA','MC']:
  s=rec['samples'][key]
  if s['min']==s['max']:extras.append(key+' constant')
  if s['nonfinite']:extras.append(key+f" nonfinite={s['nonfinite']}")
  if s['empty_visible_range']:extras.append(key+' empty in range')
 if extras:lab.SetTextSize(.024);lab.DrawLatex(.18,.60,'; '.join(extras))
 bot.cd();axis=R.TH1D(n+'_axis','',15,lo,hi);axis.SetDirectory(0);axis.SetMinimum(0);axis.SetMaximum(5);axis.GetYaxis().SetTitle('DATA / MC');axis.GetYaxis().SetTitleSize(.09);axis.GetYaxis().SetLabelSize(.08);axis.GetYaxis().SetTitleOffset(.7);axis.GetYaxis().SetNdivisions(304);axis.GetXaxis().SetTitle(cfg['title']);axis.GetXaxis().SetTitleSize(.10);axis.GetXaxis().SetLabelSize(.09);axis.Draw('AXIS')
 g=R.TGraphErrors();g.SetName(n+'_ratio');g.SetMarkerStyle(20);g.SetMarkerSize(.6);clip=0
 for i in range(1,16):
  den=mc.GetBinContent(i);num=data.GetBinContent(i)
  if den<=0:continue
  val=num/den;err=float(np.hypot(data.GetBinError(i)/den,num*mc.GetBinError(i)/den**2));k=g.GetN();g.SetPoint(k,mc.GetBinCenter(i),val);g.SetPointError(k,0,err);clip+=int(val<0 or val>5);rec['ratio_points'].append({'bin':i,'value':val,'error':err})
 g.Draw('P SAME');line=R.TLine(lo,1,hi,1);line.SetLineStyle(2);line.Draw();lab.SetTextSize(.075)
 if clip:lab.DrawLatex(.18,.82,f'{clip} points outside [0,5]')
 if not g.GetN():lab.DrawLatex(.18,.82,'No ratio: MC bins empty')
 rec['ratio_points_outside_0_5']=clip;root.cd();g.Write();c.SaveAs(str(out/'Bu/pdf'/f'{n}.pdf'));c.Close();metrics[n]=rec
root.Close();subprocess.run(['pdfunite',*[str(out/'Bu/pdf'/f'{n}.pdf') for n in names],str(out/'Bu/variables_all.pdf')],check=True)
with uproot.open(out/'metadata/histograms.root') as f:
 for n,rec in metrics.items():
  for key,s in rec['samples'].items():
   h=f[n+'_'+key];raw=f[n+'_'+key+'_raw'];assert np.allclose(raw.values(),s['raw_bins']);assert np.allclose(raw.variances(),s['sumw2_bins']);assert np.allclose(h.values(),np.array(s['raw_bins'])/s['normalization']);assert np.allclose(h.variances(),np.array(s['sumw2_bins'])/s['normalization']**2)
for p in list((out/'Bu/pdf').glob('*.pdf'))+[out/'Bu/variables_all.pdf']:
 info=subprocess.check_output(['pdfinfo',str(p)],text=True);assert int(re.search(r'Pages:\s+(\d+)',info)[1])==(len(names) if p.name=='variables_all.pdf' else 1)
report['channels']['Bu']='complete, numerical validation PASS';report['variables_Bu']=len(names);report['validation']={'ROOT_selection_sumw_principal_bins':'PASS','all_variable_flow_sumw_closure':'PASS','persisted_raw_normalized_bins_sumw2':'PASS','PDF_page_counts':'PASS','visual_review':'pending'}
for name,value in [('provenance.json',report),('variable_metrics.json',metrics),('variables.json',configs)]: (out/'metadata'/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
for p in [Path(__file__),repo/'scripts/workflows/ppref_psi2s_prompt_plots.py']:shutil.copy2(p,out/'metadata'/p.name)
(out/'README.md').write_text('# Bvar_ppRef_dataMC_2D4\n\n基础fiducial：Bpt > 7.5 GeV、abs(By) < 2.4；再加Bnorm_svpvDistance_2D > 4。无mass cut。DATA为原始含背景候选，不做拟合或减背景。MC要求pthat >= pThats、pThatreweight > 0，以pThatreweight填图，保留来源5。\n\n- [Bu合订本](Bu/variables_all.pdf) · [单变量PDF](Bu/pdf/)\n- Bs/Bd：当前新版ppRef MC未交付，仅有OLD_OLD旧文件；本次未替代使用，等待用户通知新版输入。\n\n沿用Variables_ppRef的变量顺序、范围与15 bins。各分布在可见有限范围内分别单位面积归一化，误差sqrt(sumw2)/归一化常数；不含归一化协方差。橙色MC误差带、蓝色raw DATA点及DATA/MC ratio；MC bin为零不画ratio，超出[0,5]的点数标注。范围外比例、非有限值及恒定字段逐变量记录。Bu第二径迹等原始占位字段仅作输入诊断，不推定其物理有效性。\n\n[输入和验证](metadata/provenance.json) · [逐变量计数、范围、ratio](metadata/variable_metrics.json)。源输入size/mtime与ML inventory及读前后核对，未重算源文件SHA。\n')
print('BU NUMERICAL VALIDATION PASS; Bs/Bd pending current ppRef MC',flush=True)
