from pathlib import Path
import gzip,json,hashlib,time,resource
from collections import defaultdict,Counter
ROOT=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct4-hla-allele-support');t0=time.time()
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
with gzip.open(ROOT/'common-haplotype-markers.json.gz','rt') as f:j=json.load(f)
exclude={m['sequence'] for m in json.load(open(ROOT/'nonB-genomic-mask.json'))['new_nonB_matches']}
ms=[m for m in j['markers'] if m['k']==91 and m['B_gene_exclusive'] and m['passes_transcript_paralog_mask'] and m['sequence'] not in exclude]
positions={t:{p for m in ms if m['allele']==t for rid,p in m['positions']} for t in j['targets']};shared=set.intersection(*positions.values());ms=[m for m in ms if any(p in shared for rid,p in m['positions'])];words={t:set() for t in j['targets']}
for m in ms:words[m['allele']].update([m['sequence'],rc(m['sequence'])])
def fastq(f):
 while True:
  n=f.readline()
  if not n:return
  s=f.readline().strip();plus=f.readline();q=f.readline().strip()
  assert n.startswith('@') and plus.startswith('+') and len(s)==len(q)
  yield n.split()[0][1:],s,q
co=Counter();families=defaultdict(set);names=defaultdict(list);example_support=[];thresholds={f'Q{q}_end{e}':Counter() for q in (20,25,30) for e in (0,5)};nr=0
with gzip.open(ROOT/'common-positive.R1.fastq.gz','rt') as f1,gzip.open(ROOT/'common-positive.R2.fastq.gz','rt') as f2:
 for a,b in zip(fastq(f1),fastq(f2),strict=True):
  assert a[0].removesuffix('/1')==b[0].removesuffix('/2');nr+=1;hits=[]
  for t,ww in words.items():
   for mate,(n,s,q) in enumerate([a,b]):
    for w in ww:
     start=s.find(w)
     while start>=0:
      hits.append((t,min(ord(x)-33 for x in q[start:start+91]),min(start,len(s)-start-91),mate,start));start=s.find(w,start+1)
  for q in (20,25,30):
   for e in (0,5):
    aa=sorted({t for t,qual,edge,_,_ in hits if qual>=q and edge>=e});lab=';'.join(aa) if aa else 'no_marker';thresholds[f'Q{q}_end{e}'][lab]+=1
    if q==30 and e==5 and len(aa)==1:
     t=aa[0];families[t].add(hashlib.sha256((a[1]+'|'+b[1]).encode()).hexdigest());names[t].append(hashlib.sha256(a[0].removesuffix('/1').encode()).hexdigest())
expected=json.load(open(ROOT/'common-patient-counts.json'));check=thresholds['Q30_end5']==expected['totals']['Q30_end5']['91']['matched_positions'];fcheck={t:len(v) for t,v in families.items()}=={t:expected['distinct_full_pair_sequences_matched_Q30_end5'][t]['91'] for t in j['targets']};assert check and fcheck
out={'status':'passed','method':'Separate implementation: Python gzip four-line parser, str.find scanning every marker orientation, independent quality/read-end checks; no Aho-Corasick or pysam matcher. Recounts the exhaustive full-scan positive-pair extraction, not a separate search for omitted negatives. Applies extra all-genomic non-B HLA mask.','candidate_pairs_recounted':nr,'matched_start_positions':len(shared),'threshold_counts':thresholds,'distinct_full_pair_sequences_Q30_end5':{t:len(v) for t,v in families.items()},'matches_primary_counts':check,'matches_primary_sequence_families':fcheck,'supporting_pair_name_hashes_Q30_end5':names,'elapsed_seconds':time.time()-t0,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(ROOT/'independent-recount91.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='supporting_pair_name_hashes_Q30_end5'}),flush=True)
