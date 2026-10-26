from pathlib import Path
import ROOT as R,numpy as np,uproot,json,hashlib,shutil,os
R.gROOT.SetBatch(True);R.gStyle.SetOptStat(0)
o=Path('/eos/home-l/leyao/pbpb_work/X_analysis/Analysis_CODES/plotER/results/pb18_an_score');meta_dir=o/'metadata';old=json.loads((meta_dir/'manifest.json').read_text());src=Path(old['source']);pdf=o/'Bmass_score0p06_Qmumu_5MeV.pdf';assert not pdf.exists();st=src.stat();assert st.st_size==old['source_size'] and st.st_mtime_ns==old['source_mtime_ns']
cols=['Bmass','Bpt','By','Btktkmass','BQvalue','Score']+[f'Bmu{i}{k}' for i in [1,2] for k in ['pt','eta','phi']]
with uproot.open(src) as f:a=f['ntmix'].arrays(cols,library='np')
assert all(np.isfinite(z).all() for z in a.values());vec=[];mu_mass=.10565837
for i in [1,2]:
 pt=a[f'Bmu{i}pt'].astype('float64');eta=a[f'Bmu{i}eta'].astype('float64');phi=a[f'Bmu{i}phi'].astype('float64')
 vec.append(tuple(v.astype('float32') for v in [pt*np.cos(phi),pt*np.sin(phi),pt*np.sinh(eta),np.sqrt((pt*np.cosh(eta))**2+mu_mass**2)]))
sums=[(vec[0][i]+vec[1][i]).astype('float64') for i in range(4)];mass2=sums[3]**2-sum(z*z for z in sums[:3]);assert (mass2>0).all();mm=np.sqrt(mass2).astype('float32');q=a['Bmass'].astype('float64')-mm.astype('float64')-a['Btktkmass'].astype('float64');sel=(a['Score']>.06)&(q<.13);sel_old=(a['Score']>.06)&(a['BQvalue']<.13)
# Confirm the previously matched original candidates without rereading the large legacy input.
pairs=json.loads((meta_dir/'run2_field_pairs.json').read_text());lookup={tuple(float(a[k][i]) for k in ['Bmass','Bpt','By']):i for i in range(1000)}
for row in pairs:
 i=lookup[tuple(row['flat'][k] for k in ['Bmass','Bpt','By'])];assert float(mm[i])==row['orig']['Bmumumass']
# Independent C++/TLorentzVector implementation and ROOT histogram.
R.gInterpreter.Declare('''
float pb18_restore_mumu(float pt1,float eta1,float phi1,float pt2,float eta2,float phi2){
 TLorentzVector x,y,z;x.SetPtEtaPhiM(pt1,eta1,phi1,0.10565837);y.SetPtEtaPhiM(pt2,eta2,phi2,0.10565837);
 float ax=x.Px(),ay=x.Py(),az=x.Pz(),ae=x.E(),bx=y.Px(),by=y.Py(),bz=y.Pz(),be=y.E();
 z.SetPxPyPzE(ax+bx,ay+by,az+bz,ae+be);return (float)z.M();}
''')
df=R.RDataFrame('ntmix',str(src)).Define('Bmumumass_recovered','pb18_restore_mumu(Bmu1pt,Bmu1eta,Bmu1phi,Bmu2pt,Bmu2eta,Bmu2phi)').Define('Q_mumu','double(Bmass)-double(Bmumumass_recovered)-double(Btktkmass)')
check=df.Range(1000).AsNumpy(['Bmumumass_recovered']);assert np.array_equal(check['Bmumumass_recovered'],mm[:1000])
filtered=df.Filter('Score > 0.06 && Q_mumu < 0.13');hist=filtered.Histo1D(('Bmass_Qmumu',';m_{J/#psi #pi^{+}#pi^{-}} [GeV/c^{2}];Entries / 5 MeV/c^{2}',80,3.6,4.0),'Bmass');h=hist.GetValue();edges=np.linspace(3.6,4.0,81);cnt=np.histogram(a['Bmass'][sel],edges)[0];cnt[-1]-=int(np.sum(sel&(a['Bmass']==4)))
assert int(h.GetEntries())==int(sel.sum());assert np.array_equal(cnt,[h.GetBinContent(i) for i in range(1,81)])
c=R.TCanvas('cPb18Qmumu','',1000,650);c.SetLeftMargin(.10);c.SetRightMargin(.36);c.SetTopMargin(.06);c.SetBottomMargin(.12)
h.SetLineColor(R.kBlue+1);h.SetLineWidth(2);h.SetFillColor(R.kBlue);h.SetFillStyle(3358);h.SetMinimum(0);h.SetMaximum(1.25*float(cnt.max()));h.Draw('HIST')
l=R.TLegend(.66,.49,.99,.91);l.SetBorderSize(0);l.SetFillColor(R.kWhite);l.SetTextSize(.030);l.SetHeader('#bf{PbPb 2018 Data}');l.AddEntry(h,f'N = {int(sel.sum()):,}','f')
for label in ['Existing REflated sample','Score (BDT) > 0.06','Q = m_{J/#psi#pi#pi} - m_{#mu#mu} - m_{#pi#pi}','Q < 0.13 GeV/c^{2}','m_{#mu#mu} recovered from muons']:l.AddEntry(R.nullptr,label,'')
l.Draw();c.SaveAs(str(pdf))
root=meta_dir/'Bmass_score0p06_Qmumu_5MeV.root';assert not root.exists();f=R.TFile.Open(str(root),'RECREATE');h.Write();f.Close()
recovered=meta_dir/'Qmumu_recovered.root';assert not recovered.exists()
with uproot.recreate(recovered) as f:f['recovered']={'source_entry':np.arange(len(mm),dtype='int64'),'Bmumumass_recovered':mm,'Q_mumu':q,'passes_score_and_Q':sel}
np.savetxt(meta_dir/'Bmass_bins_Qmumu_5MeV.csv',np.column_stack([edges[:-1],edges[1:],cnt,np.sqrt(cnt)]),delimiter=',',header='mass_low_GeV,mass_high_GeV,count,sqrt_count',comments='')
m={'source':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'selection':'Score > 0.06 && Q_mumu < 0.13','Q_definition':'double(Bmass) - double(Bmumumass_recovered) - double(Btktkmass)','mumu_method':'muon pt/eta/phi with m_mu=0.10565837; float32 four-vector components and summed components, double invariant mass calculation then float32 mass as original Bntuple','input_entries':len(mm),'selected_entries':int(sel.sum()),'previous_BQ_selected':int(sel_old.sum()),'common':int((sel&sel_old).sum()),'added_vs_BQ':int((sel&~sel_old).sum()),'removed_vs_BQ':int((~sel&sel_old).sum()),'mass_range':[3.6,4.0],'bin_width_GeV':.005,'bins':80,'underflow':h.GetBinContent(0),'overflow':h.GetBinContent(81),'validation':{'original_Bmumumass_exact_matches':len(pairs),'independent_ROOT_mumu_matches':1000,'ROOT_vs_numpy_all_bins':'PASS'},'scope':'Only Score and reconstructed Q cuts added to existing REflated selection; no additional AN fiducial/centrality/PVz cuts, no fit','recovered_tree':str(recovered),'recovered_tree_sha256':hashlib.sha256(recovered.read_bytes()).hexdigest(),'pdf_sha256':hashlib.sha256(pdf.read_bytes()).hexdigest()};assert m['source_sha256']==old['source_sha256'] and src.stat().st_mtime_ns==st.st_mtime_ns
(meta_dir/'Qmumu_manifest.json').write_text(json.dumps(m,indent=2)+'\n');shutil.copy2('/tmp/plot_pb18_qmumu.py',meta_dir/'plot_pb18_qmumu.py');print(json.dumps(m,indent=2))
