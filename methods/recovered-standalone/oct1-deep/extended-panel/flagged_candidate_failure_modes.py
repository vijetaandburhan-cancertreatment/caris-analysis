"""Why two Caris-flagged candidate SNPs fail a stronger local check."""
import pathlib,json,collections,statistics,pysam
B=pathlib.Path(__file__).resolve().parent;ROOT=B.parents[2];SRC=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853');out=[]
with pysam.AlignmentFile(str(SRC/'DNA_TN26-279853.bam'),'rb',index_filename=str(ROOT/'work/oct1-analysis/DNA_TN26-279853.bam.bai')) as bam:
 for gene,ch,pos,alt in [('ERBB2','chr17',39725388,'T'),('PTEN','chr10',87961098,'C')]:
  p=pos-1;ct=collections.Counter();ends=[];patterns=set()
  for r in bam.fetch(ch,p,p+1):
   if r.flag&(4|256|512|1024|2048) or r.mapping_quality<20 or r.query_qualities is None:continue
   aps={rp:q for q,rp in r.get_aligned_pairs(matches_only=True)};q=aps.get(p)
   if q is None or r.query_qualities[q]<20 or r.query_sequence[q]!=alt:continue
   ct['SNP_alt_reads']+=1;ct['reverse' if r.is_reverse else 'forward']+=1;ct['softclipped']+=int(any(op==4 for op,n in r.cigartuples));dist=min(q-r.query_alignment_start,r.query_alignment_end-1-q);ends.append(dist);patterns.add((r.reference_start,r.reference_end,r.cigarstring,r.is_reverse,r.is_read1))
   if p-8 not in aps or p+8 not in aps:ct['fails_span_8bp_flanks']+=1
   else:
    a,b=aps[p-8],aps[p+8]+1
    if min(r.query_qualities[a:b])<20:ct['fails_flank_BQ20']+=1
    else:ct['high_quality_spanning_but_other_haplotype']+=1
  out.append({'gene':gene,'pos1':pos,'counts':dict(ct),'median_distance_to_alignment_end':statistics.median(ends) if ends else None,'max_distance_to_alignment_end':max(ends) if ends else None,'distinct_alignment_patterns':len(patterns)})
(B/'flagged-candidate-failure-modes.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
