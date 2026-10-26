from pathlib import Path
import numpy as np,uproot,awkward as ak,json,hashlib
root=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872')
flatpath=root/'legacy_run2/REflated/flat_ntmix_PbPb18_MC_X3872.root'
rawpath=root/'legacy_run2/inputs/mc_prompt/prompt_x3872.root'
with uproot.open(flatpath) as file: flat=file['ntmix_X3872'].arrays(['Bmass','Bpt','By','pThatreweight'],library='np')
fields=['Bmass','Bpt','By'];dtype=np.dtype([(name,'<f4') for name in fields])
def signatures(values):
 out=np.empty(len(values['Bmass']),dtype=dtype)
 for key in fields:out[key]=values[key]
 return out
keys=signatures(flat);order=np.argsort(keys);sortedkeys=keys[order];matches=np.zeros(len(keys),dtype=np.int32);weight_matches=np.zeros(len(keys),dtype=np.int32);source_out=np.full(len(keys),-1,dtype=np.int32);pthat_out=np.full(len(keys),np.nan);event_out=np.full(len(keys),-1,dtype=np.int64)
with uproot.open(rawpath) as file:
 event=file['hiEvtAnalyzer/HiTree'].arrays(['sample','pthat','pthatweight'],library='np');offset=0
 for batch in file['Bfinder/ntmix'].iterate(fields,step_size=100000,library='ak'):
  flatraw={key:ak.to_numpy(ak.flatten(batch[key])) for key in fields};rawkeys=signatures(flatraw)
  lo=np.searchsorted(sortedkeys,rawkeys,side='left');hi=np.searchsorted(sortedkeys,rawkeys,side='right');hit=np.flatnonzero(hi>lo)
  event_indices=ak.to_numpy(ak.flatten(ak.broadcast_arrays(batch.Bmass,ak.Array(np.arange(offset,offset+len(batch))))[1]))
  for idx in hit:
   targets=order[lo[idx]:hi[idx]];origin=event_indices[idx];matches[targets]+=1;weight_matches[targets]+=(flat['pThatreweight'][targets]==event['pthatweight'][origin]);source_out[targets]=int(event['sample'][origin]);pthat_out[targets]=float(event['pthat'][origin]);event_out[targets]=origin
  offset+=len(batch);print('events read',offset,flush=True)
report={'raw':str(rawpath),'flat':str(flatpath),'signature':fields,'flat_n':len(keys),'no_match':int((matches==0).sum()),'ambiguous':int((matches>1).sum()),'exactly_one_match':int((matches==1).sum()),'no_equal_weight_match':int((weight_matches==0).sum()),'source5_candidates':int((source_out==5).sum()),'zero_weight_candidates':int((flat['pThatreweight']==0).sum()),'raw_event_count':offset}
assert (matches>=1).all() and (weight_matches>=1).all(),report
assert (matches==1).all();np.savez_compressed('/tmp/pb18_reco_source_map.npz',sample=source_out,pthat=pthat_out,source_event=event_out);report['status']='PASS';Path('/tmp/run2_reco_inheritance.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
