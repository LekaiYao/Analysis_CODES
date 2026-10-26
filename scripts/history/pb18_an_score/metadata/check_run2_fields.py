import uproot,numpy as np,json
from pathlib import Path
base=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872/legacy_run2')
with uproot.open(base/'REflated/flat_ntmix_PbPb18_DATA.root') as f:a=f['ntmix'].arrays(['Bmass','Bpt','By','CentBin','Bujmass','BQvalue','Btktkmass','PVz'],entry_stop=1000,library='np')
keys={tuple(float(a[k][i]) for k in ['Bmass','Bpt','By']):i for i in range(1000)};out=[]
with uproot.open(base/'inputs/data/nominal_data_flatten.root') as f:
 for c in f['ntmix_flatten'].iterate(['Bmass','Bpt','By','hiBin','Bujmass','Bmumumass','Btktkmass','PVz'],step_size=500000,library='np'):
  for j in np.flatnonzero(np.isin(c['Bmass'],a['Bmass'])):
   key=tuple(float(c[k][j]) for k in ['Bmass','Bpt','By'])
   if key in keys:
    i=keys.pop(key);out.append({'flat':{k:float(a[k][i]) for k in a},'orig':{k:float(c[k][j]) for k in c}})
  if not keys:break
Path('/tmp/run2_field_pairs.json').write_text(json.dumps(out,indent=2))
print('matched',len(out));print(json.dumps(out[:6],indent=2))
print('cent_equal',sum(x['flat']['CentBin']==x['orig']['hiBin'] for x in out),'cent_half',sum(x['flat']['CentBin']==x['orig']['hiBin']/2 for x in out))
