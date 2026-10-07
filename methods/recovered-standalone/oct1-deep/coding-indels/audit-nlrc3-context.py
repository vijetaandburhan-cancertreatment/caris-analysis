from pathlib import Path
from collections import Counter,defaultdict
import json,pysam
P=Path(__file__).resolve().parent;ROOT=P.parents[2];RAW=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
w=json.loads((P/'reference-windows/chr16-3563534-3563695.json').read_text());ref=w['dna'].upper();pos=3563614
allout=[]
for kind in ('DNA','RNA'):
 b=pysam.AlignmentFile(str(RAW/f'{kind}_TN26-279853.bam'),'rb',index_filename=str(ROOT/f'work/oct1-analysis/{kind}_TN26-279853.bam.bai'))
 haplotypes=defaultdict(set);detail={};other=defaultdict(set)
 for r in b.fetch('chr16',pos-80,pos+81):
  if r.flag&(4|256|512|1024|2048) or r.mapping_quality<20 or (r.has_tag('NH') and r.get_tag('NH')>1):continue
  name=r.query_name;events=[];rp=r.reference_start;qp=0;covered=[]
  for op,n in r.cigartuples:
   if op in (0,7,8):
    for j in range(max(0,pos-30-rp),min(n,pos+32-rp)):
     if r.query_qualities[qp+j]>=20:
      covered.append(rp+j);observed=r.query_sequence[qp+j];expected=ref[rp+j-w['start']]
      if observed!=expected:events.append(('SNV',rp+j,expected,observed))
    qp+=n;rp+=n
   elif op==1:
    if pos-30<=rp<pos+32 and min(r.query_qualities[qp:qp+n])>=20:events.append(('I',rp,r.query_sequence[qp:qp+n]))
    qp+=n
   elif op==2:
    if pos-30<=rp<pos+32:events.append(('D',rp,n))
    rp+=n
   elif op==3:rp+=n
   elif op==4:qp+=n
  if not (set(range(pos-12,pos+13)) | {3563596})<=set(covered):continue
  key=tuple(events);haplotypes[key].add(name)
  if kind=='RNA':detail.setdefault(str(key),[]).append({'name':name,'cigar':r.cigarstring,'start0':r.reference_start,'end0':r.reference_end,'reverse':r.is_reverse})
 b.close()
 rows=[{'events':[list(a) for a in ev],'paired_names':len(names)} for ev,names in haplotypes.items()]
 rows.sort(key=lambda r:r['paired_names'],reverse=True);allout.append({'kind':kind,'haplotypes_observed_with_25bp_core_and_neighbor_SNV_covered':rows,'RNA_details':detail})
(P/'audit-nlrc3-context.json').write_text(json.dumps(allout,indent=2)+'\n')
for r in allout:print(r['kind'],r['haplotypes_observed_with_25bp_core_and_neighbor_SNV_covered'][:8])
