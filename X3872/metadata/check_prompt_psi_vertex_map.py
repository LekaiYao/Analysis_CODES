from pathlib import Path
import uproot,numpy as np,json
root=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872');raw=root/'legacy_run2/inputs/mc_prompt/prompt_psi2s.root';flat=root/'legacy_run2/REflated/flat_ntmix_PbPb18_MC_PSI2S.root'
fields=['PVx','PVy','PVz'];dtype=np.dtype([(key,'<f4') for key in fields])
def pack(values):
 output=np.empty(len(values['PVz']),dtype=dtype)
 for key in fields:output[key]=values[key]
 return output
with uproot.open(raw) as file:
 print({key:file['Bfinder/ntmix'][key].typename for key in fields},flush=True)
 events=file['Bfinder/ntmix'].arrays(fields,library='np');weights=file['hiEvtAnalyzer/HiTree'].arrays(['sample','pthat','pthatweight'],library='np')
with uproot.open(flat) as file:candidates=file['ntmix_PSI2S'].arrays(fields+['pThatreweight'],library='np')
keys=pack(events);candidatekeys=pack(candidates);order=np.argsort(keys);sortedkeys=keys[order];lower=np.searchsorted(sortedkeys,candidatekeys,'left');upper=np.searchsorted(sortedkeys,candidatekeys,'right');counts=upper-lower
print('counts',np.unique(counts,return_counts=True),flush=True)
assert (counts==1).all()
origins=order[lower];assert np.array_equal(candidates['pThatreweight'],weights['pthatweight'][origins]);np.savez_compressed('/tmp/prompt_psi_vertex_map.npz',pThats=weights['sample'][origins].astype(np.int32),pthat=weights['pthat'][origins],source_event=origins,old_weight=weights['pthatweight'][origins]);print('Unique PVxyz event mapping and old weights PASS',flush=True)
