from pathlib import Path
import pysam,json,bisect,collections,sys,time,hashlib

def panel(path):
 d=collections.defaultdict(set)
 for ln in Path(path).read_text().splitlines():
  ch,s,e=ln.split('\t');assert int(e)==int(s)+1;d[ch].add(int(s))
 return {ch:sorted(s) for ch,s in d.items()}
def matches(rd,refs,sites,names):
 named=rd.query_name in names
 overlap=False
 if not rd.is_unmapped:
  ch=refs[rd.reference_id];p=sites.get(ch,())
  i=bisect.bisect_left(p,rd.reference_start)
  if i<len(p) and p[i]<rd.reference_end:
   for s,e in rd.get_blocks():
    j=bisect.bisect_left(p,s)
    if j<len(p) and p[j]<e:overlap=True;break
 return named,overlap

def capture(bam,bed,qnames,out,summary,maxbytes=256*1024**2):
 sites=panel(bed);names=set(Path(qnames).read_text().splitlines());t=time.time();counts=collections.Counter();seen=set();flag=collections.Counter();mapped=collections.Counter()
 with pysam.AlignmentFile(str(bam),'rb') as inp,pysam.AlignmentFile(str(out),'wb',template=inp,threads=2) as dest:
  refs=inp.references
  for rd in inp.fetch(until_eof=True):
   counts['input_records']+=1;named,overlap=matches(rd,refs,sites,names)
   if named or overlap:
    dest.write(rd);counts['captured_records']+=1;counts['named_records']+=named;counts['cohort_overlap_records']+=overlap;counts['both_records']+=named and overlap
    if named:
     seen.add(rd.query_name);flag[str(rd.flag)]+=1;mapped['unmapped' if rd.is_unmapped else refs[rd.reference_id]]+=1
   if counts['input_records']%1000000==0:
    info={'counts':dict(counts),'exception_names_seen':len(seen),'elapsed_seconds':time.time()-t,'output_bytes':Path(out).stat().st_size};Path(str(summary)+'.progress').write_text(json.dumps(info,indent=2)+'\n')
    if info['output_bytes']>maxbytes:raise RuntimeError('Capture size bound exceeded; input never fully saved')
  assert Path(out).stat().st_size<=maxbytes
 missing=sorted(names-seen)
 info={'status':'CAPTURE_STREAM_EOF','counts':dict(counts),'frozen_exception_names':len(names),'exception_names_seen':len(seen),'missing_exception_names':missing,'named_record_flags':dict(flag),'named_record_locations':dict(mapped),'elapsed_seconds':time.time()-t,'output_bytes':Path(out).stat().st_size,'selection':'OR of fixed SNP M/=/X base overlap and exact frozen QNAME. No MAPQ/BQ/flag filtering at capture; all emitted records for named reads across all locations retained. Not all possible alignments.','bam_header':inp.header.to_dict() if False else 'Preserved in BAM'}
 Path(summary).write_text(json.dumps(info,indent=2)+'\n');print(json.dumps({k:info[k] for k in ['status','counts','exception_names_seen','output_bytes']},indent=2),flush=True)

if __name__=='__main__':capture(*sys.argv[1:6])
