#!/usr/bin/env python3
import json,pathlib,pysam,hashlib,collections
from Bio.Seq import Seq
H=pathlib.Path.home();B=H/'.local/share/codex/caris-analysis';O=pathlib.Path(__file__).parent;D=B/'oct5-fusion-priority-review'
audit=json.loads((B/'oct4-fusion-read-audit/patient/junction-read-audit.json').read_text());out={}
f=pysam.FastaFile(str(B/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'));bam=pysam.AlignmentFile(str(B/'TN26-279853/DNA_TN26-279853.bam'),'rb',index_filename=str(pathlib.Path.cwd()/'work/oct1-analysis/DNA_TN26-279853.bam.bai'))
rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def key(a):return(a.query_name,a.get_tag('RG')if a.has_tag('RG')else'')
for ci in ['078bc453229c','e2a8e3625460']:
 r=next(r for r in audit['rows']if r['candidate_id']==ci);chrom,s1=r['breakpoint1'].split(':');p1=int(s1);p2=int(r['breakpoint2'].split(':')[1]);seen=set();reads=[]
 intervals=[(p2-501,p1+500)]if p1-p2<1000 else[(p1-501,p1+500),(p2-501,p2+500)]
 for s,e in intervals:
  for a in bam.fetch(chrom,s,e):
   kk=(a.query_name,a.flag,a.reference_start,a.cigarstring)
   if kk in seen:continue
   seen.add(kk)
   if not a.flag&3844 and a.mapping_quality>=20 and a.query_sequence and a.query_qualities:reads.append(a)
 l,rr=r['fusion_transcript'].split('|');l=l.upper();rr=rr.upper();assert set(l+rr)<=set('ACGT')
 z={'gene':r['gene1'],'longer_exact_junction_markers':{}}
 for arm in(25,35,50,65,75):
  if len(l)<arm or len(rr)<arm:continue
  pat=l[-arm:]+rr[:arm]
  for q in(20,30):
   pairs=collections.defaultdict(list);seqfamilies=set();alnfamilies=set();n_fwd=n_rev=0
   for a in reads:
    s=a.query_sequence;i=s.find(pat)
    # BAM sequences are reference-oriented; these self events are plus-strand.
    if i<0 or min(a.query_qualities[i:i+len(pat)])<q:continue
    pairs[key(a)].append(a);seqfamilies.add(min(s,rc(s)));alnfamilies.add((a.reference_start,a.reference_end,a.cigarstring,a.is_read1,a.is_reverse));n_fwd+=not a.is_reverse;n_rev+=a.is_reverse
   z['longer_exact_junction_markers'][f'{arm}+{arm}_Q{q}']={'fragments':len(pairs),'distinct_whole_read_sequences_orientation_normalized':len(seqfamilies),'distinct_alignment_descriptors':len(alnfamilies),'forward_alignments':n_fwd,'reverse_alignments':n_rev,'both_mates':sum(len({a.is_read1 for a in vv})==2 for vv in pairs.values())}
 if ci=='e2a8e3625460':
  pos=41190343;insert='GCTGTGGGTCCAGCTGCTGCCAGCCTA';z['insertion_after_pos1']=pos;z['inserted_sequence']=insert
  # exact genomic-coordinate flanks, full interval from first to last aligned query base, no N operation allowed.
  calls={}
  for flank in(10,20):
   s=pos-flank;e=pos+flank;ref=f.fetch(chrom,s,e).upper();alt=ref[:flank]+insert+ref[flank:];names=collections.defaultdict(set);readsamples=collections.Counter()
   for a in reads:
    pairs={rp:qp for qp,rp in a.get_aligned_pairs(matches_only=True)if rp is not None};indices=[pairs.get(rp)for rp in range(s,pos)]+[pairs.get(rp)for rp in range(pos,e)]
    if None in indices:continue
    left=indices[:flank];right=indices[flank:]
    if any(y!=x+1 for x,y in zip(left,left[1:]))or any(y!=x+1 for x,y in zip(right,right[1:])):continue
    start=left[0];end=right[-1]+1
    if end<=start or min(a.query_qualities[start:end])<20:continue
    ss=a.query_sequence[start:end];cl='alt'if ss==alt else'ref'if ss==ref else'other';names[key(a)].add(cl);readsamples[cl]+=1
   cc=collections.Counter(next(iter(v))if len(v)==1 else'conflicting_mates'for v in names.values())
   calls[f'exact{flank}bp_flanks_all_BQ20_MAPQ20']={'fragments':dict(cc),'read_alignments':dict(readsamples),'reference_haplotype':ref,'alternate_haplotype':alt}
  z['coordinate_anchored_haplotype_counts']=calls
  s,e=41190320,41190410;ref=f.fetch(chrom,s,e).upper();dna=ref[:pos-s]+insert+ref[pos-s:]
  ll=list(f.fetch(chrom,s,41190374).upper());ll[41190370-1-s]='A';rna=''.join(ll)+f.fetch(chrom,41190348-1,e).upper();assert dna==rna
  z['caller_ITD_plus_mismatch_equals_DNA_27bp_insertion']=True
  cds=f.fetch(chrom,41189887-1,41190636).upper();altcds=cds[:pos-(41189887-1)]+insert+cds[pos-(41189887-1):];wt=str(Seq(cds).translate());mut=str(Seq(altcds).translate());prefix=next(i for i,(a,b)in enumerate(zip(wt,mut))if a!=b);suffix=0
  while suffix<min(len(wt)-prefix,len(mut)-prefix) and wt[-suffix-1]==mut[-suffix-1]:suffix+=1
  z['GENCODE37_ENST00000398470_1_translation']={'reference_coding_bp':len(cds),'alternate_coding_bp':len(altcds),'reference_aa':len(wt),'alternate_aa':len(mut),'shared_prefix_aa':prefix,'shared_suffix_aa':suffix,'reference_changed_segment':wt[prefix:len(wt)-suffix],'alternate_changed_segment':mut[prefix:len(mut)-suffix],'reference_local_aa':wt[prefix-12:len(wt)-suffix+12],'alternate_local_aa':mut[prefix-12:len(mut)-suffix+12],'internal_stop':('*'in mut)}
 if ci=='078bc453229c':z['GENCODE37_ENST00000415495_5_coding_end']='stop codon chr19:41807218-41807220; donor41807383 is downstream in3primeUTR; no translated junction predicted'
 out[ci]=z
out['limits']=['Local DNA evidence corroborates sequence but does not establish somatic status; no matched normal.','CEACAM3 longer markers corroborate the joined sequence, but repeat/paralog alternatives and exact genomic structural interpretation need independent validation.','KRTAP9-1 coordinate-anchored count excludes clipped/alternative representations and is not a VAF.','Sequencing duplicates not always flagged; unique read names/whole-read sequence/alignment descriptors are not independent original molecules.','KRTAP translated sequence is a selected reference isoform prediction, not protein detection or clinical significance.']
for folder in(O,D):(folder/'bounded-DNA-details.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out))
