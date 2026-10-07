from pathlib import Path
import sys,json,gzip,hashlib,time,resource,itertools
from collections import defaultdict,Counter
ROOT=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct4-hla-allele-support');sys.path.insert(0,str(ROOT/'deps'));import ahocorasick,pysam
RAW=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853');paths=[RAW/f'RNA_TN26-279853_S25.R{x}.fastq.gz' for x in (1,2)];t0=time.time()
with gzip.open(ROOT/'common-haplotype-markers.json.gz','rt') as f:j=json.load(f)
markers=[m for m in j['markers'] if m['B_gene_exclusive'] and m['passes_transcript_paralog_mask']];targets=list(j['targets']);KS=[31,51,71,91,111,131]
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
pos={t:{k:{p for m in markers if m['allele']==t and m['k']==k for rid,p in m['positions']} for k in KS} for t in targets}
matched={k:pos[targets[0]][k]&pos[targets[1]][k] for k in KS}
for m in markers:m['matched_position']=any(p in matched[m['k']] for rid,p in m['positions'])
a=ahocorasick.Automaton();idx=defaultdict(list)
for i,m in enumerate(markers):
 for w in {m['sequence'],rc(m['sequence'])}:idx[w].append(i)
for w,ids in idx.items():a.add_word(w,(len(w),ids))
a.make_automaton()
# Public equal-start/same-position controls and naive matcher parity before any patient scan.
control=[]
for t,r in j['targets'].items():
 s=r['sequence']
 for seq in [s,rc(s)]:
  fast=sorted((p-k+1,i) for p,(k,ids) in a.iter(seq) for i in ids)
  slow=sorted((p,i) for w,ids in idx.items() for p in range(len(seq)-len(w)+1) if seq[p:p+len(w)]==w for i in ids)
  assert fast==slow
  control.append({'allele':t,'reverse':seq!=s,'naive_parity':True})
opportunities=[]
for k in KS:
 for L in (75,100,125,150):
  for fragment in (200,250,300):
   passed=[]
   for t,r in j['targets'].items():
    s=r['sequence'];total=hit=0
    for start in range(len(s)-fragment+1):
     spans=[(start+5,start+L-5),(start+fragment-L+5,start+fragment-5)];total+=1
     hit+=any(any(lo<=p and p+k<=hi for p in matched[k]) for lo,hi in spans)
    passed.append(hit);opportunities.append({'allele':t,'k':k,'read_length':L,'fragment_length':fragment,'end_distance':5,'total_starts':total,'marker_positive_starts':hit})
   assert passed[0]==passed[1],(k,L,fragment,passed)
(ROOT/'common-public-controls.json').write_text(json.dumps({'exact_match_controls':control,'equal_start_opportunities':opportunities,'matched_start_counts':{str(k):len(v) for k,v in matched.items()},'passed':True},indent=2))
levels=[(q,e) for q in (20,25,30) for e in (0,5)];totals={f'Q{q}_end{e}':{str(k):{'all_gene_specific':Counter(),'matched_positions':Counter()} for k in KS} for q,e in levels};joint={f'Q{q}_end{e}':Counter() for q,e in levels};counts=Counter();mcount={f'Q{q}_end{e}':Counter() for q,e in levels};families=defaultdict(set);combined_groups=Counter();positives=[gzip.open(ROOT/f'common-positive.R{x}.fastq.gz','wt',compresslevel=5) for x in (1,2)];evidence=gzip.open(ROOT/'common-positive-evidence.jsonl.gz','wt',compresslevel=5)
def norm(s):return s[:-2] if s.endswith(('/1','/2')) else s
with pysam.FastxFile(str(paths[0])) as f1,pysam.FastxFile(str(paths[1])) as f2:
 for r1,r2 in itertools.zip_longest(f1,f2):
  if r1 is None or r2 is None:raise RuntimeError('unequal EOF')
  if norm(r1.name)!=norm(r2.name):raise RuntimeError('mate name mismatch')
  counts['pairs']+=1;records=[r1,r2];hits=[]
  for mate,r in enumerate(records):
   if len(r.sequence)!=len(r.quality):raise RuntimeError('quality length mismatch')
   for end,(k,ids) in a.iter(r.sequence):
    start=end-k+1;q=min(ord(c)-33 for c in r.quality[start:end+1]);edge=min(start,len(r.sequence)-end-1)
    hits.extend((i,mate,start,q,edge) for i in ids)
  if hits:
   counts['any_matching_pair']+=1
   for h,r in zip(positives,records):h.write(f'@{r.name}'+(' '+r.comment if r.comment else '')+'\n'+r.sequence+'\n+\n'+r.quality+'\n')
   evidence.write(json.dumps({'name_hash_sha256':hashlib.sha256(norm(r1.name).encode()).hexdigest(),'hits':hits})+'\n')
   for q,e in levels:
    key=f'Q{q}_end{e}';ids={i for i,m,p,b,d in hits if b>=q and d>=e};mcount[key].update(ids)
    for k in KS:
     these={i for i in ids if markers[i]['k']==k}
     for mode in ('all_gene_specific','matched_positions'):
      aa=sorted({markers[i]['allele'] for i in these if mode=='all_gene_specific' or markers[i]['matched_position']})
      label=';'.join(aa) if aa else 'no_marker';totals[key][str(k)][mode][label]+=1
      if key=='Q30_end5' and mode=='matched_positions' and len(aa)==1:
       families[(aa[0],k)].add(hashlib.sha256((r1.sequence+'|'+r2.sequence).encode()).hexdigest())
    byallele={t:{i for i in ids if markers[i]['allele']==t} for t in targets}
    present=[t for t in targets if byallele[t]]
    if len(present)==1:
     t=present[0];groups=set.intersection(*(set(markers[i]['compatible_two_field']) for i in byallele[t]))
     joint[key][t+';'+('unique_two_field_intersection' if groups=={t} else 'multiple_B_alleles')]+=1
     if key=='Q30_end5':combined_groups[(t,tuple(sorted(groups)))]+=1
    elif len(present)>1:joint[key]['cross_allele_conflict']+=1
  if counts['pairs']%1000000==0:
   prog={'stage':'counting_common_haplotype_markers','pairs':counts['pairs'],'positive_pairs':counts['any_matching_pair'],'elapsed_seconds':time.time()-t0,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss};(ROOT/'common-patient-progress.json').write_text(json.dumps(prog,indent=2));print(json.dumps(prog),flush=True)
for h in positives:h.close()
evidence.close();assert counts['pairs']==23209264
out={'status':'passed','counts':dict(counts),'inputs':[str(p) for p in paths],'totals':totals,'joint_competing_group_intersections':joint,'joint_groups_Q30_end5':[{'allele':t,'groups':list(g),'pairs':n} for (t,g),n in combined_groups.most_common()],'distinct_full_pair_sequences_matched_Q30_end5':{t:{str(k):len(families[(t,k)]) for k in KS} for t in targets},'markers':markers,'marker_pair_counts':{k:dict(v) for k,v in mcount.items()},'elapsed_seconds':time.time()-t0,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'normal_EOF_and_mate_name_and_seq_quality_lengths':True,'method':'Exact31/51/71/91/111/131mer matches to common B*50:01:01/B*52:01:01 coding representatives, either orientation; mask matches to other HLA genes in all available IMGT3.46 CDS and non-HLA-B GENCODE37 transcripts; count query-name pairs at minimum whole-marker Q20/25/30,0/5nt end distance. Matched windows retain only shared coding-start coordinates admissible for both allele sequences. This is conditional expression evidence given prior HLA typing.','limits':['No individual common-reference marker through131nt is specific to a single two-field allele across every competing HLA-B sequence.','Joint within-fragment intersections exclude known conflicting sequence combinations but cannot exclude unknown or incompletely sequenced alleles.','No calibrated expression ratio, tumor-cell identity, allele retention/LOH, copy number or protein/peptide presentation measurement.','Read-pair and distinct full sequence counts are not independently tagged molecules.']}
(ROOT/'common-patient-counts.json').write_text(json.dumps(out,indent=2));print(json.dumps({'status':out['status'],'counts':out['counts'],'Q30_end5':totals['Q30_end5'],'joint_Q30_end5':joint['Q30_end5'],'distinct_pairs':out['distinct_full_pair_sequences_matched_Q30_end5'],'elapsed_seconds':out['elapsed_seconds'],'peak_RSS_bytes':out['peak_RSS_bytes']}),flush=True)
