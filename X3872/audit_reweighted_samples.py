from pathlib import Path
from decimal import Decimal,localcontext
import argparse,csv,hashlib,json
import numpy as np
import uproot
import ROOT as R

parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True);args=parser.parse_args();root=args.root.resolve();metadata=root/'metadata'
manifest=json.loads((metadata/'mc_manifest.json').read_text());assert manifest['status']=='PASS' and manifest['source5_policy']=='zero' and len(manifest['files'])==8
R.gROOT.SetBatch(True)

def digest(path):
 checksum=hashlib.sha256()
 with path.open('rb') as stream:
  for chunk in iter(lambda:stream.read(8388608),b''):checksum.update(chunk)
 return checksum.hexdigest()

csv_path=metadata/'Henrique_mc_normalization_4855f01.csv';assert digest(csv_path)==manifest['2023_normalization_table']['sha256'];table=list(csv.DictReader(line for line in csv_path.read_text().splitlines() if not line.startswith('#')))
records=[];boundaries=[10,15,30,50,float('inf')]
for record in manifest['files']:
 path=Path(record['output']);path=path if path.is_absolute() else root.parent/path;source=Path(record['source']);assert path.exists();assert digest(path)==record['output_sha256'];stat=source.stat();assert (stat.st_size,stat.st_mtime_ns)==(record['source_size'],record['source_mtime_ns']);assert digest(source)==record['source_sha256']
 with localcontext() as context:
  context.prec=45
  if record['year']=='2023':
   rows=sorted([row for row in table if (row['system'],row['particle'],row['promptness'])==('PbPb23',record['particle'],record['kind'])],key=lambda row:int(row['pthat']));denominator=Decimal(0);expected=[]
   for row in rows[1:]:
    denominator+=Decimal(row['n_gen'])/(Decimal(row['xsec_pb'])*Decimal(row['filter_eff']));expected.append(float(1/denominator))
  else:
   with uproot.open(source) as input_file:
    reco_name=next(name for name in input_file.keys(cycle=False) if name!='ntGen');original_weights=np.unique(input_file[reco_name]['pThatreweight'].array(library='np'));positive=sorted([float(weight) for weight in original_weights if weight>0],reverse=True);assert len(positive)==5
   low=Decimal.from_float(positive[0]);expected=[float(1/(1/Decimal.from_float(weight)-1/low)) for weight in positive[1:]]
 root_file=R.TFile.Open(str(path));file_result={'path':str(path),'year':record['year'],'particle':record['particle'],'kind':record['kind'],'positive_interval_weights':expected,'trees':{}}
 with uproot.open(path) as incoming,uproot.open(source) as original:
  assert incoming.keys(cycle=False)==original.keys(cycle=False)
  for name in incoming.keys(cycle=False):
   tree=incoming[name];columns=['pthat','pThats','pThatreweight']+(['Bpt','By','BQvalue'] if name!='ntGen' else []);values=tree.arrays(columns,library='np');actual=values['pthat'];component=values['pThats'];weight=values['pThatreweight'];assert tree.num_entries==original[name].num_entries
   zero=(component==5)|(actual<component);assert np.array_equal(weight==0,zero);assert np.isfinite(weight).all()
   intervals=[]
   for index,(lower,upper) in enumerate(zip(boundaries[:-1],boundaries[1:])):
    mask=(actual>=lower)&(actual<upper)&(~zero);assert np.allclose(weight[mask],expected[index],rtol=3e-13,atol=0)
    intervals.append({'lower':lower,'upper':str(upper),'count':int(mask.sum()),'weight':expected[index]})
   root_tree=root_file.Get(name);assert root_tree.GetEntries('pThatreweight == 0')==int(zero.sum());assert root_tree.GetEntries('pThats==5 && pThatreweight!=0')==0;assert root_tree.GetEntries('pthat<pThats && pThatreweight!=0')==0
   result={'entries':tree.num_entries,'positive':int((~zero).sum()),'zero':int(zero.sum()),'source5':int((component==5).sum()),'below_own_threshold':int((actual<component).sum()),'valid_intervals':intervals}
   if name!='ntGen':
    fid=(values['Bpt']>15)&(values['Bpt']<50)&(abs(values['By'])<1.6)&(values['BQvalue']<.15);selected=weight[fid];fid_positive=fid&(~zero);result['fiducial']={'expression':'Bpt>15 && Bpt<50 && abs(By)<1.6 && BQvalue<0.15','candidates':int(fid.sum()),'positive_candidates':int(fid_positive.sum()),'zero_candidates':int((fid&zero).sum()),'sumw':float(selected.sum()),'sumw2':float(selected@selected),'neff':float(selected.sum()**2/(selected@selected))}
    assert root_tree.GetEntries('Bpt>15 && Bpt<50 && abs(By)<1.6 && BQvalue<0.15 && pThatreweight>0')==int(fid_positive.sum())
   file_result['trees'][name]=result
 root_file.Close();file_result['status']='PASS';records.append(file_result);print('Independent audit PASS',record['year'],record['particle'],record['kind'],flush=True)
for record in json.loads((metadata/'data_copies.json').read_text()):
 assert digest(Path(record['output']))==record['sha256'];assert digest(Path(record['source']))==record['sha256']
assert not list(root.rglob('*.partial'))
output={'status':'PASS','MC_files':8,'DATA_byte_copies':2,'rows_preserved_all_trees':True,'source5_and_below_threshold_zero_exact':True,'independent_decimal_weights':'PASS','ROOT_zero_and_fiducial_counts':'PASS','data_hashes':'PASS','records':records,'audit_script_sha256':digest(Path(__file__))}
(metadata/'independent_acceptance.json').write_text(json.dumps(output,indent=2)+'\n');print('ALL INDEPENDENT CHECKS PASS',flush=True)
