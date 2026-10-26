import uproot,numpy as np,json
from pathlib import Path
p='/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2/REflated/flat_ntmix_PbPb18_DATA.root'
with uproot.open(p) as f:a=f['ntmix'].arrays(['Bmass','Bpt','By','Bmu1pt','Bmu2pt','Bmu1eta','Bmu2eta','Bmu1phi','Bmu2phi'],entry_stop=1000,library='np')
a={k:v.astype('float64') for k,v in a.items()};m=.10565837
v=[]
for i in [1,2]:
 pt=a[f'Bmu{i}pt'];eta=a[f'Bmu{i}eta'];phi=a[f'Bmu{i}phi'];v.append((pt*np.cos(phi),pt*np.sin(phi),pt*np.sinh(eta),np.sqrt((pt*np.cosh(eta))**2+m*m)))
mm=np.sqrt(np.maximum(0,(v[0][3]+v[1][3])**2-sum((v[0][i]+v[1][i])**2 for i in range(3))))
vf=[tuple(z.astype('float32') for z in vv) for vv in v]
sums=[(vf[0][i]+vf[1][i]).astype('float64') for i in range(4)]
mm_float=np.sqrt(np.maximum(0,sums[3]**2-sum(z*z for z in sums[:3]))).astype('float32')
d_float=[]
lookup={tuple(float(a[k][i]) for k in ['Bmass','Bpt','By']):i for i in range(1000)};d=[];examples=[]
for row in json.load(open('/tmp/run2_field_pairs.json')):
 o=row['orig'];i=lookup[tuple(row['flat'][k] for k in ['Bmass','Bpt','By'])];d.append(mm[i]-o['Bmumumass']);d_float.append(float(mm_float[i])-o['Bmumumass'])
 if len(examples)<3:examples.append({'mumu_reconstructed':mm[i],'original_Bmumumass':o['Bmumumass'],'Bujmass':row['flat']['Bujmass']})
r={'matched':len(d),'muon_mass_GeV':m,'max_abs_delta_GeV':float(np.max(np.abs(d))),'original_float_intermediate_max_delta_GeV':float(np.max(np.abs(d_float))),'original_float_intermediate_exact_matches':int(np.sum(np.array(d_float)==0)),'rms_delta_GeV':float(np.sqrt(np.mean(np.array(d)**2))),'examples':examples};assert len(d)==1000;r['within_2e_6_GeV']=bool(np.max(np.abs(d))<2e-6);Path('/tmp/run2_mumu_reconstruction.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
