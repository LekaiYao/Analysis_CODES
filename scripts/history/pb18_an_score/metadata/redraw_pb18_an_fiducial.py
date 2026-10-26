from pathlib import Path
import ROOT as R,numpy as np,uproot,json,hashlib,shutil,os,tempfile
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0)
o=Path('/eos/home-l/leyao/pbpb_work/X_analysis/Analysis_CODES/plotER/results/pb18_an_score');meta=o/'metadata';stage=Path(tempfile.mkdtemp(prefix='pb18_an_fid_'));(stage/'metadata').mkdir()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
m=json.loads((meta/'manifest.json').read_text());qm=json.loads((meta/'Qmumu_manifest.json').read_text());src=Path(m['source']);assert sha(src)==m['source_sha256'];assert sha(meta/'Qmumu_recovered.root')==qm['recovered_tree_sha256']
with uproot.open(src) as f:a=f['ntmix'].arrays(['Bmass','Bpt','By','PVz','Score','BQvalue'],library='np')
with uproot.open(meta/'Qmumu_recovered.root') as f:q=f['recovered'].arrays(library='np')
assert np.array_equal(q['source_entry'],np.arange(len(a['Bmass'])))
fid=(a['Bpt']>15)&(a['Bpt']<50)&(abs(a['By'])<1.6)&(abs(a['PVz'])<15)
selection='Bpt > 15 && Bpt < 50 && abs(By) < 1.6 && abs(PVz) < 15'
R.gInterpreter.Declare('''float restore_pb18_mm(float pt1,float eta1,float phi1,float pt2,float eta2,float phi2){TLorentzVector x,y,z;x.SetPtEtaPhiM(pt1,eta1,phi1,0.10565837);y.SetPtEtaPhiM(pt2,eta2,phi2,0.10565837);float ax=x.Px(),ay=x.Py(),az=x.Pz(),ae=x.E(),bx=y.Px(),by=y.Py(),bz=y.Pz(),be=y.E();z.SetPxPyPzE(ax+bx,ay+by,az+bz,ae+be);return (float)z.M();}''')
df=R.RDataFrame('ntmix',str(src)).Define('mm_recovered','restore_pb18_mm(Bmu1pt,Bmu1eta,Bmu1phi,Bmu2pt,Bmu2eta,Bmu2phi)').Define('Q_mumu','double(Bmass)-double(mm_recovered)-double(Btktkmass)')
results={};arrays={}
for label,width,name,hname,csvname in [('BQ',5,'Bmass_score0p06','Bmass_score0p06','Bmass_bins'),('BQ',10,'Bmass_score0p06_10MeV','Bmass_score0p06_10MeV','Bmass_bins_10MeV'),('Q',5,'Bmass_score0p06_Qmumu_5MeV','Bmass_Qmumu','Bmass_bins_Qmumu_5MeV')]:
 cut='Score > 0.06 && '+('BQvalue < 0.13' if label=='BQ' else 'Q_mumu < 0.13')+' && '+selection
 mask=(a['Score']>.06)&((a['BQvalue']<.13) if label=='BQ' else (q['Q_mumu']<.13))&fid
 n=int(mask.sum());assert n==(10528 if label=='BQ' else 10424)
 bins=400//width;edges=np.linspace(3.6,4.,bins+1);count=np.histogram(a['Bmass'][mask],edges)[0];count[-1]-=int(np.sum(mask&(a['Bmass']==4.)))
 ptr=df.Filter(cut).Histo1D((hname,';m_{J/#psi #pi^{+}#pi^{-}} [GeV/c^{2}];Entries / '+str(width)+' MeV/c^{2}',bins,3.6,4.),'Bmass');h=ptr.GetValue();assert int(h.GetEntries())==n;assert np.array_equal(count,[h.GetBinContent(i) for i in range(1,bins+1)]);assert np.allclose(count,[h.GetBinError(i)**2 for i in range(1,bins+1)]);assert h.GetBinContent(0)==h.GetBinContent(bins+1)==0
 c=R.TCanvas('c_'+name,'',1000,650);c.SetLeftMargin(.10);c.SetRightMargin(.36);c.SetTopMargin(.06);c.SetBottomMargin(.12)
 h.SetLineColor(R.kBlue+1);h.SetLineWidth(2);h.SetFillColor(R.kBlue);h.SetFillStyle(3358);h.SetMinimum(0);h.SetMaximum(1.25*float(count.max()));h.Draw('HIST')
 l=R.TLegend(.66,.35,.995,.91);l.SetBorderSize(0);l.SetFillColor(R.kWhite);l.SetTextSize(.029);l.SetHeader('#bf{PbPb 2018 Data}');l.AddEntry(h,f'N = {n:,}','f')
 lines=['Score (BDT) > 0.06']
 if label=='BQ':lines+=['BQvalue < 0.13 GeV/c^{2}','Fixed J/#psi mass in BQvalue']
 else:lines+=['Q = m_{J/#psi#pi#pi} - m_{#mu#mu} - m_{#pi#pi}','Q < 0.13 GeV/c^{2}','m_{#mu#mu} recovered from muons']
 lines+=['15 < Bpt < 50 GeV/c','|By| < 1.6','|PVz| < 15 cm']
 for text in lines:l.AddEntry(R.nullptr,text,'')
 l.Draw();c.SaveAs(str(stage/(name+'.pdf')))
 f=R.TFile.Open(str(stage/'metadata'/(name+'.root')),'RECREATE');h.Write();f.Close()
 np.savetxt(stage/'metadata'/(csvname+'.csv'),np.column_stack([edges[:-1],edges[1:],count,np.sqrt(count)]),delimiter=',',header='mass_low_GeV,mass_high_GeV,count,sqrt_count',comments='')
 results[name]={'selection':cut,'selected_entries':n,'bins':bins,'bin_width_GeV':width/1000,'ROOT_numpy_bins_and_counts':'PASS','errors':'sqrt(count)','pdf_sha256':sha(stage/(name+'.pdf'))};arrays[name]=count
assert np.array_equal(arrays['Bmass_score0p06'].reshape(40,2).sum(axis=1),arrays['Bmass_score0p06_10MeV'])
m.update(selection=results['Bmass_score0p06']['selection'],selected_entries=10528,histogram_entries=10528,scope='BDT and BQ cuts plus AN Bpt, By, PVz; centrality not reapplied',additional_an_cuts_not_applied={'centrality_percent':[0,90]},pdf_sha256=results['Bmass_score0p06']['pdf_sha256'],plot_script='metadata/redraw_pb18_an_fiducial.py',plot_script_sha256=sha(Path(__file__)),user_decision='2026-09-29: add AN Bpt/By/PVz and overwrite Q/BQ plots',visual_review='pending')
qm.update(selection=results['Bmass_score0p06_Qmumu_5MeV']['selection'],selected_entries=10424,previous_BQ_selected=10528,scope='BDT and recovered Q plus AN Bpt, By, PVz; centrality not reapplied',pdf_sha256=results['Bmass_score0p06_Qmumu_5MeV']['pdf_sha256'],plot_script='metadata/redraw_pb18_an_fiducial.py',plot_script_sha256=sha(Path(__file__)))
bq=(a['Score']>.06)&(a['BQvalue']<.13)&fid;qq=(a['Score']>.06)&(q['Q_mumu']<.13)&fid
qm.update(common=int((bq&qq).sum()),added_vs_BQ=int((qq&~bq).sum()),removed_vs_BQ=int((bq&~qq).sum()));qm['validation']['visual_review']='pending';qm['recovered_tree_note']='Cached passes_score_and_Q only records Score/Q; current plot additionally applies Bpt/By/PVz by source_entry. Restored masses unchanged.'
rebin={'source_histogram':str(meta/'Bmass_score0p06.root'),'source_sha256':sha(stage/'metadata/Bmass_score0p06.root'),'selection':m['selection'],'bin_width_GeV':.01,'bins':40,'range':[3.6,4.0],'entries':10528,'verification':'ROOT vs NumPy all bins and errors PASS; 40 bins equal adjacent 5 MeV bin sums','pdf_sha256':results['Bmass_score0p06_10MeV']['pdf_sha256']}
for name,data in [('manifest.json',m),('Qmumu_manifest.json',qm),('rebin_10MeV.json',rebin),('fiducial_redraw.json',results)]: (stage/'metadata'/name).write_text(json.dumps(data,indent=2)+'\n')
(stage/'README.md').write_text('''# PbPb18：AN score、Q与运动学选择后的质量谱

**当前三张图已加入15 < Bpt < 50、abs(By) < 1.6、abs(PVz) < 15 cm。**

| Q选择（均要求Score > 0.06及上述三项cut） | 候选数 | 图 |
|---|---:|---|
| BQvalue < 0.13 | 10528 | [5 MeV/bin](Bmass_score0p06.pdf) · [10 MeV/bin](Bmass_score0p06_10MeV.pdf) |
| Bmass − m(μμ) − Btktkmass < 0.13 | 10424 | [5 MeV/bin](Bmass_score0p06_Qmumu_5MeV.pdf) |

BQvalue使用固定J/ψ质量3.096916 GeV；原文Q使用逐候选m(μμ)。m(μμ)由flat的muon四动量恢复，先前与1000个原始Bmumumass逐值核验一致。原恢复缓存不变，其中passes_score_and_Q只含Score/Q；本次绘图额外施加三项cut。

输入仍为Henrique的2018 REflated DATA，单位计数；质量范围3.6–4.0 GeV，沿用原图布局和样式。Score与原始BDT抽查10000条一致，AN2019_021_v10表3给出BDT>0.06。中心度编码尚未确认，未追加中心度cut；不宣称完整复现AN全部选择。未拟合。

ROOT与NumPy独立选择计数、全部质量bin及误差通过，10 MeV每bin等于相邻两个5 MeV bin之和。当前可重复入口为metadata/redraw_pb18_an_fiducial.py；其它旧绘图脚本保留历史来源含义。

[本次验收](metadata/fiducial_redraw.json) · [BQ来源及选择](metadata/manifest.json) · [Q来源及选择](metadata/Qmumu_manifest.json) · [10 MeV核验](metadata/rebin_10MeV.json)
''')
shutil.copy2(__file__,stage/'metadata/redraw_pb18_an_fiducial.py')
# Validate and back up every existing artifact before replacement.
files=[p for p in stage.rglob('*') if p.is_file()];backup=meta/'originals/before_fiducial';assert not backup.exists();backup.mkdir();old={}
for p in files:
 rel=p.relative_to(stage);dest=o/rel
 if dest.exists():old[str(rel)]=sha(dest);target=backup/rel;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(dest,target)
for p in files:
 rel=p.relative_to(stage);dest=o/rel
 if str(rel) in old:assert sha(dest)==old[str(rel)]
 else:assert not dest.exists()
 tmp=dest.with_name('.'+dest.name+'.fid_tmp');assert not tmp.exists();shutil.copy2(p,tmp)
 if dest.exists():shutil.copymode(dest,tmp)
 os.replace(tmp,dest)
shutil.rmtree(stage)
print(json.dumps(results,indent=2));print('PASS: all three plots replaced, originals recoverable')
