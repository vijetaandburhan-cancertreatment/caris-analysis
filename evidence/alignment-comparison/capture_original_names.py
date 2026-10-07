from pathlib import Path
import json,pysam,time,hashlib,collections,shutil,sys
R=Path(__file__).resolve().parent;O=R/'original-complete-names';O.mkdir(exist_ok=True)
assert json.loads((R/'patient-full-pass/run.json').read_text())['status']=='COMPLETE'
names=set((R/'frozen-exception-qnames.txt').read_text().splitlines());assert len(names)==3426
source=R.parent/'TN26-279853/RNA_TN26-279853.bam';dest=O/'original-3426-names.bam';assert not dest.exists()
c=collections.Counter();seen=set();t=time.time();minfree=shutil.disk_usage(O).free
with pysam.AlignmentFile(str(source),'rb') as inp,pysam.AlignmentFile(str(dest),'wb',template=inp,threads=2) as out:
 for rd in inp.fetch(until_eof=True):
  c['input_records']+=1
  if rd.query_name in names:
   out.write(rd);seen.add(rd.query_name);c['captured_records']+=1
  if c['input_records']%1000000==0:
   free=shutil.disk_usage(O).free;minfree=min(minfree,free)
   assert free>=3.5*1024**3,'Disk floor'
   assert dest.stat().st_size<64*1024**2,'Output bound'
   assert time.time()-t<900,'Time bound'
pysam.quickcheck(str(dest));assert seen==names
with dest.open('rb') as f:h=hashlib.file_digest(f,'sha256').hexdigest()
r={'status':'COMPLETE','source':str(source),'output':str(dest),'output_SHA256':h,'source_selection':'Exact frozen3426QNAMEs, all records across original full BAM. No location, flag, MAPQ, BQ restriction.','counts':dict(c),'names_seen':len(seen),'missing_names':sorted(names-seen),'elapsed_seconds':time.time()-t,'output_bytes':dest.stat().st_size,'min_free_disk_bytes':minfree};(O/'summary.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
