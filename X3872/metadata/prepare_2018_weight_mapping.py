from pathlib import Path
import argparse,json,hashlib
import numpy as np,uproot,awkward as ak
parser=argparse.ArgumentParser();parser.add_argument('particle',choices=['x3872','psi2s']);parser.add_argument('kind',choices=['prompt','nonprompt']);args=parser.parse_args()
root=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872');label='X3872' if args.particle=='x3872' else 'PSI2S';suffix='' if args.kind=='prompt' else '_nonPrompt'
raw=root/'legacy_run2/inputs'/('mc_'+args.kind)/(args.kind+'_'+args.particle+'.root');flat=root/'legacy_run2/REflated'/('flat_ntmix_PbPb18_MC_'+label+suffix+'.root');output=Path('/tmp/pthat_rebuild_maps')/(args.kind+'_'+args.particle);output.mkdir(parents=True,exist_ok=True)
assert not (output/'acceptance.json').exists()
stat_raw=raw.stat();stat_flat=flat.stat()
with uproot.open(flat) as file:
 reco=file['ntmix_'+label].arrays(['Bmass','Bpt','By','pThatreweight'],library='np');gen_columns=[name for name in file['ntGen'].keys() if name!='pThatreweight'];gen=file['ntGen'].arrays(gen_columns,library='np')
fields=['Bmass','Bpt','By'];dtype=np.dtype([(name,'<f4') for name in fields])
def signatures(values):
 result=np.empty(len(values['Bmass']),dtype=dtype)
 for name in fields:result[name]=values[name]
 return result
keys=signatures(reco);order=np.argsort(keys);ordered=keys[order];matches=np.zeros(len(keys),dtype=np.int32);weight_matches=np.zeros(len(keys),dtype=np.int32);sources=np.full(len(keys),-1,dtype=np.int32);actual=np.full(len(keys),np.nan,dtype=np.float32);event_index=np.full(len(keys),-1,dtype=np.int64)
with uproot.open(raw) as file:
 event=file['hiEvtAnalyzer/HiTree'].arrays(['sample','pthat','pthatweight'],library='np');event_count=len(event['sample']);assert event_count==file['Bfinder/ntmix'].num_entries==file['Bfinder/ntGen'].num_entries
 offset=0
 if args.particle=='x3872' and args.kind=='prompt':
  known=np.load('/eos/home-l/leyao/pbpb_work/X_analysis/Analysis_CODES/plotER/results/prompt_x_pthat5/metadata/pb18_reco_source_map.npz');sources[:]=known['sample'];actual[:]=known['pthat'];event_index[:]=known['source_event'];assert np.array_equal(sources,event['sample'][event_index]);assert np.array_equal(actual,event['pthat'][event_index]);assert np.array_equal(reco['pThatreweight'],event['pthatweight'][event_index]);matches[:]=1;weight_matches[:]=1
 else:
  for batch in file['Bfinder/ntmix'].iterate(fields,step_size=100000,library='ak'):
   flattened={name:ak.to_numpy(ak.flatten(batch[name])) for name in fields};raw_keys=signatures(flattened);lower=np.searchsorted(ordered,raw_keys,side='left');upper=np.searchsorted(ordered,raw_keys,side='right');hits=np.flatnonzero(upper>lower)
   origins=ak.to_numpy(ak.flatten(ak.broadcast_arrays(batch.Bmass,ak.Array(np.arange(offset,offset+len(batch))))[1]))
   for index in hits:
    targets=order[lower[index]:upper[index]];origin=origins[index];matches[targets]+=1;weight_matches[targets]+=(reco['pThatreweight'][targets]==event['pthatweight'][origin]);sources[targets]=int(event['sample'][origin]);actual[targets]=event['pthat'][origin];event_index[targets]=origin
   offset+=len(batch);print(args.kind,args.particle,'RECO events',offset,flush=True)
 assert (matches==1).all() and (weight_matches==1).all(),(int((matches==0).sum()),int((matches>1).sum()),int((weight_matches!=1).sum()))
 np.savez_compressed(output/'reco.npz',pThats=sources,pthat=actual,source_event=event_index,old_weight=event['pthatweight'][event_index])
 offset=0;candidate_offset=0;parts={name:[] for name in ['pThats','pthat','source_event','old_weight']}
 for batch in file['Bfinder/ntGen'].iterate(gen_columns+['GisSignal'],step_size=100000,library='ak'):
  selected=(batch.GisSignal==7)&(batch.Gpt>=5)&(abs(batch.Gy)<=2.4);size=int(ak.sum(selected));end=candidate_offset+size
  for name in gen_columns:
   values=ak.to_numpy(ak.flatten(batch[name][selected]));assert np.array_equal(values,gen[name][candidate_offset:end],equal_nan=True),(name,candidate_offset,end)
  origins=ak.to_numpy(ak.flatten(ak.broadcast_arrays(batch.Gpt,ak.Array(np.arange(offset,offset+len(batch))))[1][selected]))
  parts['source_event'].append(origins)
  for name,column in [('pThats','sample'),('pthat','pthat'),('old_weight','pthatweight')]:parts[name].append(event[column][origins])
  candidate_offset=end;offset+=len(batch)
 assert candidate_offset==len(gen['Gpt']);assert offset==event_count
 np.savez_compressed(output/'gen.npz',**{name:np.concatenate(values) for name,values in parts.items()})
 result={'status':'PASS','particle':args.particle,'kind':args.kind,'raw':{'path':str(raw),'size':stat_raw.st_size,'mtime_ns':stat_raw.st_mtime_ns,'uuid':str(file.file.uuid)},'flat':{'path':str(flat),'size':stat_flat.st_size,'mtime_ns':stat_flat.st_mtime_ns},'reco_candidates':len(keys),'GEN_candidates':candidate_offset,'reco_unique_match':True,'reco_old_weights_exact':True,'GEN_all_original_columns_exact':gen_columns,'GEN_selection':'GisSignal==7 && Gpt>=5 && abs(Gy)<=2.4','files':{path.name:hashlib.sha256(path.read_bytes()).hexdigest() for path in [output/'reco.npz',output/'gen.npz']}}
assert raw.stat().st_mtime_ns==stat_raw.st_mtime_ns and flat.stat().st_mtime_ns==stat_flat.st_mtime_ns
(output/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
