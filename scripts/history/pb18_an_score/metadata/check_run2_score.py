from pathlib import Path
import uproot,numpy as np,json
flat=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_DATA.root');orig=flat.parent.parent/'inputs/data/nominal_data_flatten.root'
with uproot.open(flat) as f:a=f['ntmix'].arrays(['Bmass','Bpt','By','Score','CentBin','BQvalue','Bujmass','Btktkmass'],library='np')
idx=np.linspace(0,len(a['Bmass'])-1,10000,dtype=int);keys={tuple(float(a[k][i]) for k in ['Bmass','Bpt','By']):int(i) for i in idx};counts={n:0 for n in ['BDT','BDTG','BDTD','BDTF']};seen=set();cent=[];q=[];examples=[]
with uproot.open(orig) as f:
 for c in f['ntmix_flatten'].iterate(['Bmass','Bpt','By','BDT','BDTG','BDTD','BDTF','hiBin','Bmumumass','Bujmass','Btktkmass'],step_size=500000,library='np'):
  # A sorted mass lookup limits Python work to exact mass matches.
  samplem=np.array([k[0] for k in keys]);possible=np.flatnonzero(np.isin(c['Bmass'],samplem))
  for j in possible:
   key=tuple(float(c[k][j]) for k in ['Bmass','Bpt','By'])
   if key not in keys:continue
   i=keys[key]
   if i in seen:continue
   seen.add(i)
   for n in counts:counts[n]+=int(a['Score'][i]==c[n][j])
   cent.append(float(a['CentBin'][i]-c['hiBin'][j]/2))
   q.append([float(a['BQvalue'][i]-(c['Bmass'][j]-c['Bmumumass'][j]-c['Btktkmass'][j])),float(a['BQvalue'][i]-(c['Bmass'][j]-3.096916-c['Btktkmass'][j]))])
   if len(examples)<5:examples.append({'Bmass':key[0],'Bpt':key[1],'Score':float(a['Score'][i]),'BDT':float(c['BDT'][j])})
r={'requested':len(keys),'matched':len(seen),'score_exact_matches':counts,'centbin_minus_hibin_half_unique':np.unique(cent).tolist(),'q_max_abs_difference_mumu_vs_fixed':np.max(np.abs(q),axis=0).tolist() if q else None,'examples':examples,'flat':str(flat),'original':str(orig)}
Path('/tmp/run2_score_identity.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
