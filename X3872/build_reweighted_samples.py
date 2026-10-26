from pathlib import Path
import argparse,csv,hashlib,json,os,subprocess,datetime
import numpy as np
import uproot
import ROOT as R

parser=argparse.ArgumentParser()
parser.add_argument('--source5-policy',required=True,choices=['zero','drop'])
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--maps',type=Path,default=Path(__file__).resolve().parent/'metadata/run2_source_maps')
parser.add_argument('--writer',type=Path,default=Path(__file__).resolve().parent/'ReweightFlatTree.C')
args=parser.parse_args()
repo=Path(__file__).resolve().parent.parent;source=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872');metadata=args.output/'metadata';metadata.mkdir(parents=True,exist_ok=True)
assert not (metadata/'mc_manifest.json').exists()
R.gROOT.SetBatch(True);R.gInterpreter.Declare(args.writer.read_text())
commit='4855f0156052a15096e04c07acb466e658940179'
csv_text=subprocess.check_output(['git','show',commit+':flatER/config/mc_normalization.csv'],cwd=repo,text=True)
rows=list(csv.DictReader(line for line in csv_text.splitlines() if not line.startswith('#')))
raw_audit=json.loads((repo/'docs/reference/upstream_review_evidence/2026-09-29/pthat_layers/pthat_layers_audit.json').read_text())
thresholds=np.array([5,10,15,30,50],dtype=np.int32)

def digest(path):
 checksum=hashlib.sha256()
 with path.open('rb') as stream:
  for chunk in iter(lambda:stream.read(8388608),b''):checksum.update(chunk)
 return checksum.hexdigest()

def new_weights(actual,component,luminosities):
 valid=(component!=5)&(actual>=component)
 assert np.isin(component,thresholds).all() and np.isfinite(actual).all()
 intervals=np.searchsorted(thresholds,actual,side='right')-1
 cumulative=np.cumsum(np.r_[0.,luminosities[1:]])
 assert ((intervals[valid]>=1)&(intervals[valid]<5)).all()
 weights=np.zeros(len(actual),dtype=np.float64);weights[valid]=1/cumulative[intervals[valid]]
 assert np.isfinite(weights).all() and (weights>=0).all()
 return weights

fixture_actual=np.array([9.9,10,14.9,15,29.9,30,49.9,50,100,20,20],dtype=np.float32)
fixture_source=np.array([10,10,15,15,30,30,50,50,50,5,10],dtype=np.int32)
fixture=new_weights(fixture_actual,fixture_source,np.array([9.,1.,2.,3.,4.]))
assert np.allclose(fixture,[0,1,0,1/3,0,1/6,0,1/10,1/10,0,1/3])
manifest={'status':'BUILDING','date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source5_policy':args.source5_policy,'below_own_source_threshold':'weight=0 when actual pthat < source threshold; equality is retained','bin_definition':'[10,15), [15,30), [30,50), [50,infinity)','2023_normalization_table':{'repository':'https://github.com/Legoinha/Analysis_CODES','commit':commit,'path':'flatER/config/mc_normalization.csv','sha256':hashlib.sha256(csv_text.encode()).hexdigest()},'2018_assumption':'User-approved inverse cumulative luminosity formula inferred from existing positive interval weights; original overall scale retained; absolute production recipe not independently established.','files':[]}
(metadata/'Henrique_mc_normalization_4855f01.csv').write_text(csv_text)
for year,sub in [('2023','PbPb23'),('2018','legacy_run2/REflated')]:
 for particle,label in [('x3872','X3872'),('psi2s','PSI2S')]:
  for kind,suffix in [('prompt',''),('nonprompt','_nonPrompt')]:
   filename=f'flat_ntmix_PbPb{year[2:]}_MC_{label}{suffix}.root';src=source/sub/filename;dst=args.output/sub/filename;dst.parent.mkdir(parents=True,exist_ok=True);temporary=dst.with_suffix('.root.partial');assert not dst.exists() and not temporary.exists();before=src.stat()
   if year=='2023':
    group=sorted([row for row in rows if (row['system'],row['particle'],row['promptness'])==('PbPb23',particle,kind)],key=lambda row:int(row['pthat']));assert [int(row['pthat']) for row in group]==thresholds.tolist();luminosities=np.array([int(row['n_gen'])/(float(row['xsec_pb'])*float(row['filter_eff'])) for row in group]);normalization={'method':'Henrique historical table, omit source5 from cumulative sum','rows':group,'source_L':luminosities.tolist()}
   else:
    original=raw_audit['raw2018'][kind+'_'+particle+'.root'];old=[]
    for interval in original['intervals'][1:]:
     values=[weight for weight in interval['weights'] if weight>0];assert len(values)==1;old.append(values[0])
    inverse=1/np.array(old);luminosities=np.diff(np.r_[0.,inverse]);assert (luminosities>0).all();assert np.allclose(1/np.cumsum(luminosities),old,rtol=1e-14)
    normalization={'method':'User-approved inversion of original positive interval weights','original_weights':old,'inferred_source_L':luminosities.tolist()}
    mapping_dir=args.maps/(kind+'_'+particle);acceptance=json.loads((mapping_dir/'acceptance.json').read_text());assert acceptance['status']=='PASS';assert acceptance['flat']['size']==before.st_size and acceptance['flat']['mtime_ns']==before.st_mtime_ns
    raw_stat=Path(acceptance['raw']['path']).stat();assert raw_stat.st_size==acceptance['raw']['size'] and raw_stat.st_mtime_ns==acceptance['raw']['mtime_ns'];normalization['mapping_acceptance']=acceptance
   file_record={'year':year,'particle':particle,'kind':kind,'source':str(src),'source_size':before.st_size,'source_mtime_ns':before.st_mtime_ns,'source_sha256':digest(src),'output':str(dst),'normalization':normalization,'trees':{}}
   root_source=R.TFile.Open(str(src));root_target=R.TFile.Open(str(temporary),'RECREATE','',root_source.GetCompressionSettings());tree_data={}
   with uproot.open(src) as incoming:
    for name in incoming.keys(cycle=False):
     tree=incoming[name];assert tree.classname=='TTree'
     if year=='2023':
      identity=tree.arrays(['pthat','pThats'],library='np');actual=identity['pthat'];component=identity['pThats']
     else:
      mapfile=mapping_dir/('gen.npz' if name=='ntGen' else 'reco.npz');assert digest(mapfile)==acceptance['files'][mapfile.name]
      with np.load(mapfile) as saved:actual=saved['pthat'].astype(np.float32);component=saved['pThats'].astype(np.int32)
     assert len(actual)==tree.num_entries
     weights=new_weights(actual,component,luminosities);keep=(component!=5) if args.source5_policy=='drop' else np.ones(len(weights),dtype=bool);tree_data[name]=(weights,actual,component,keep)
     R.WriteReweightedTree(root_source.Get(name),root_target,R.std.vector('double')(weights),R.std.vector('float')(actual),R.std.vector('int')(component),year=='2018',args.source5_policy=='drop')
     file_record['trees'][name]={'input_entries':len(weights),'output_entries':int(keep.sum()),'removed_source5_entries':int((~keep).sum()),'source5_input_entries':int((component==5).sum()),'below_own_threshold_input_entries':int((actual<component).sum()),'zero_weight_output_entries':int(((weights==0)&keep).sum()),'sumw':float(weights[keep].sum()),'sumw2':float(weights[keep]@weights[keep]),'weights_by_source':{str(int(value)):np.unique(weights[keep&(component==value)]).tolist() for value in np.unique(component[keep])}}
   root_target.Close();root_source.Close()
   with uproot.open(src) as incoming,uproot.open(temporary) as outgoing:
    assert incoming.keys(cycle=False)==outgoing.keys(cycle=False)
    for name in incoming.keys(cycle=False):
     original=incoming[name];written=outgoing[name];weights,actual,component,keep=tree_data[name];expected_types=original.typenames().copy()
     if year=='2018':expected_types.update(pthat='float',pThats='int32_t')
     assert written.classname=='TTree' and written.title==original.title and written.typenames()==expected_types
     original_arrays=original.arrays(library='np');written_arrays=written.arrays(library='np')
     for branch,values in original_arrays.items():
      if branch=='pThatreweight':assert np.array_equal(written_arrays[branch],weights[keep])
      else:assert np.array_equal(written_arrays[branch],values[keep],equal_nan=True),(filename,name,branch)
     if year=='2018':assert np.array_equal(written_arrays['pthat'],actual[keep]) and np.array_equal(written_arrays['pThats'],component[keep])
     assert (written_arrays['pThatreweight'][written_arrays['pthat']<written_arrays['pThats']]==0).all()
     if args.source5_policy=='zero':assert (written_arrays['pThatreweight'][written_arrays['pThats']==5]==0).all()
     else:assert (written_arrays['pThats']!=5).all()
     file_record['trees'][name]['unchanged_columns_exact']='PASS';file_record['trees'][name]['schema_and_new_weights_exact']='PASS'
   after=src.stat();assert before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns
   file_record['output_sha256']=digest(temporary);os.replace(temporary,dst);file_record['status']='PASS';manifest['files'].append(file_record);(metadata/'mc_manifest.partial.json').write_text(json.dumps(manifest,indent=2)+'\n');print('MC PASS',year,particle,kind,flush=True)
manifest['status']='PASS';manifest['script_sha256']=digest(Path(__file__));manifest['writer_sha256']=digest(args.writer);manifest['boundary_fixture']='PASS';(metadata/'mc_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');(metadata/'mc_manifest.partial.json').unlink();print('ALL MC PASS',flush=True)
