from pathlib import Path
import hashlib,json,os,shutil
root=Path('/eos/home-l/leyao/pbpb_work/X_analysis/Analysis_CODES/X3872');source=Path('/eos/user/h/hmarques/RUN3_Data_MC_sharing/X3872');metadata=root/'metadata';metadata.mkdir(parents=True,exist_ok=True)
assert not (metadata/'data_copies.json').exists()
records=[]
for sub,name in [('PbPb23','flat_ntmix_PbPb23_DATA.root'),('legacy_run2/REflated','flat_ntmix_PbPb18_DATA.root')]:
 src=source/sub/name;dst=root/sub/name;dst.parent.mkdir(parents=True,exist_ok=True);assert not dst.exists();temporary=dst.with_suffix('.root.partial');assert not temporary.exists();before=src.stat();checksum=hashlib.sha256()
 with src.open('rb') as incoming,temporary.open('xb') as outgoing:
  for chunk in iter(lambda:incoming.read(8388608),b''):outgoing.write(chunk);checksum.update(chunk)
 after=src.stat();assert before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns;verify=hashlib.sha256()
 with temporary.open('rb') as incoming:
  for chunk in iter(lambda:incoming.read(8388608),b''):verify.update(chunk)
 assert verify.hexdigest()==checksum.hexdigest();shutil.copystat(src,temporary);os.replace(temporary,dst)
 records.append({'source':str(src),'output':str(dst),'bytes':before.st_size,'source_mtime_ns':before.st_mtime_ns,'sha256':checksum.hexdigest(),'byte_exact_copy':True});print('DATA copy PASS',dst,flush=True)
(metadata/'data_copies.json').write_text(json.dumps(records,indent=2)+'\n')
