"""Indexed, reference-aware recount of candidate anchors; no network."""
from pathlib import Path
import csv,json,collections,math,time
import pysam
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
rows=list(csv.DictReader((OUT/'candidate-motif-context.tsv').open(),delimiter='\t'))
a=json.loads((OUT/'target-annotation.json').read_text());prov=json.loads((OUT/'reference-motif-provenance.json').read_text());bygene=collections.defaultdict(list)
for r in rows:bygene[r['gene']].append(r)
BAM='/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853/RNA_TN26-279853.bam';INDEX=str(ROOT/'work/oct1-analysis/RNA_TN26-279853.bam.bai')
results=[];sample_reads={}
for gene,group in bygene.items():
 g=a['genes'][gene];ref=json.load(open(prov['references'][gene]['path'])); dna=ref['dna'].upper(); base=ref['start'];candidates={(int(r['intron_start0']),int(r['intron_end0'])):r for r in group};data={key:collections.defaultdict(set) for key in candidates};examples=collections.defaultdict(list)
 with pysam.AlignmentFile(BAM,'rb',index_filename=INDEX) as bam:
  for read in bam.fetch(g['chrom'],g['start0'],g['end0']):
   if read.flag&(4|256|512|1024|2048) or read.mapping_quality<20 or read.has_tag('NH') and read.get_tag('NH')!=1:continue
   c=read.cigartuples
   if not c or not any(op==3 for op,n in c):continue
   pos=read.reference_start;q=0;qs=read.query_sequence;qual=read.query_qualities
   for i,(op,n) in enumerate(c):
    key=(pos,pos+n)
    if op==3 and key in candidates and 0<i<len(c)-1 and c[i-1][0] in (0,7,8) and c[i+1][0] in (0,7,8) and min(c[i-1][1],c[i+1][1])>=12 and qual is not None and min(qual[q-12:q+12])>=20:
     d=data[key];frag=(read.get_tag('RG') if read.has_tag('RG') else '')+'|'+read.query_name; sig=(read.reference_start,read.reference_end,read.is_reverse)
     d['base_fragments'].add(frag);d['base_signatures'].add(sig)
     lref=dna[pos-base-12:pos-base];rref=dna[pos+n-base:pos+n-base+12]
     left=qs[q-12:q];right=qs[q:q+12];perfect12=left==lref and right==rref
     if perfect12:d['perfect12_fragments'].add(frag);d['perfect12_signatures'].add(sig)
     perfect20=False
     if min(c[i-1][1],c[i+1][1])>=20 and min(qual[q-20:q+20])>=20:
      perfect20=qs[q-20:q]==dna[pos-base-20:pos-base] and qs[q:q+20]==dna[pos+n-base:pos+n-base+20]
      if perfect20:d['perfect20_fragments'].add(frag);d['perfect20_signatures'].add(sig)
     if len(examples[key])<15:
      examples[key].append({'query_name':read.query_name,'start0':read.reference_start,'end0':read.reference_end,'cigar':read.cigarstring,'flag':read.flag,'NH':read.get_tag('NH') if read.has_tag('NH') else None,'MAPQ':read.mapping_quality,'nM':read.get_tag('nM') if read.has_tag('nM') else None,'left_query12':left,'right_query12':right,'left_ref12':lref,'right_ref12':rref,'perfect12':perfect12,'perfect20':perfect20,'left_anchor':c[i-1][1],'right_anchor':c[i+1][1]})
    if op in (0,2,3,7,8):pos+=n
    if op in (0,1,4,7,8):q+=n
 for key,r in candidates.items():
  d=data[key];assert len(d['base_fragments'])==int(r['fragments']),(gene,key,len(d['base_fragments']),r['fragments'])
  r={**r,**{k:len(d[k]) for k in ['base_fragments','base_signatures','perfect12_fragments','perfect12_signatures','perfect20_fragments','perfect20_signatures']}}
  r['reference_aware_research_pass']=int(r['major_or_minor_canonical_motif']=='1' and r['perfect12_fragments']>=5 and r['perfect12_signatures']>=3 and r['perfect20_fragments']>=3)
  results.append(r);sample_reads[f'{gene}:{g["chrom"]}:{key[0]}-{key[1]}']={'row':r,'examples':examples[key]}
with (OUT/'candidate-reference-audit.tsv').open('w') as fh:
 w=csv.DictWriter(fh,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(sorted(results,key=lambda r:r['perfect12_fragments'],reverse=True))
(OUT/'candidate-reference-read-evidence.json').write_text(json.dumps(sample_reads,indent=2)+'\n')
summary={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'input_candidates':len(rows),'base_counts_independently_reproduced':True,'reference_aware_research_pass':sum(r['reference_aware_research_pass'] for r in results),'noncanonical_reference_motif':sum(r['major_or_minor_canonical_motif']=='0' for r in results),'method':'Independent indexed fetch by gene reproduces original qualifying RG+name counts. Adds exact public-hg38 reference match on both 12bp anchors and 20bp subset. Research pass adds reference-major/minor canonical motif and original count/diversity thresholds; not tumor-specificity or clinical validation. True nearby substitutions can reduce reference-perfect counts.'}
(OUT/'candidate-reference-audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
