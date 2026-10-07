#!/usr/bin/env python3
import pathlib,json,pysam,hashlib,collections
from Bio import Align
P=pathlib.Path(__file__).parent;B=pathlib.Path.home()/'.local/share/codex/caris-analysis'
def rc(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
a=Align.PairwiseAligner(mode='global',match_score=2,mismatch_score=-4)
for gap in ('insertion','deletion'):
 setattr(a,f'open_internal_{gap}_score',-6);setattr(a,f'extend_internal_{gap}_score',-1)
for side in ('left','right'):
 setattr(a,f'open_{side}_deletion_score',0);setattr(a,f'extend_{side}_deletion_score',0);setattr(a,f'open_{side}_insertion_score',-2);setattr(a,f'extend_{side}_insertion_score',-2)
fa=pysam.FastaFile(str(B/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'))
models=json.loads((P.parent/'proposed-models.json').read_text());model=models[0]['sequence']
names={r['name_sha256'] for r in json.loads((P/'softclip-SA-mate-details.json').read_text()) if any(g['SA_other_locus'] or g['mate_other_locus'] for g in r['geometry'])}
rows=[];group=collections.defaultdict(list)
with pysam.AlignmentFile(str(P/'bounded-DNA-records.bam'),'rb')as bam:
 for x in bam.fetch(until_eof=True):
  n=hashlib.sha256(x.query_name.encode()).hexdigest()
  if n not in names:continue
  row={'name_sha256':n,'flag':x.flag,'MAPQ':x.mapping_quality,'chrom':x.reference_name,'start0':x.reference_start,'CIGAR':x.cigarstring,'sequence':x.query_sequence,'qualities':list(x.query_qualities),'tags':dict(x.tags),'reverse':x.is_reverse,'read1':x.is_read1,'duplicate':x.is_duplicate,'secondary':x.is_secondary,'supplementary':x.is_supplementary,'mate_chrom':x.next_reference_name,'mate_start0':x.next_reference_start,'TLEN':x.template_length}
  rows.append(row)
  if not(x.is_secondary or x.is_supplementary):group[x.query_sequence].append(row)
refs=[{'id':'proposed_identical750nt','sequence':model,'junction_index0':250}]
for label,chrom,start in [('CCND1','chr11',69650841),('HEG1','chr3',125013197)]:
 refs.append({'id':label+'_contiguous750nt','sequence':fa.fetch(chrom,start,start+750).upper(),'chrom':chrom,'start0':start,'end0':start+750})
results=[];export=[]
for qi,(q,members)in enumerate(group.items(),1):
 for ori in ('+','-'):
  query=q if ori=='+'else rc(q)
  for ref in refs:
   al=a.align(ref['sequence'],query)[0];co=al.coordinates.tolist();mm=[];ins=[];dels=[];paired=[]
   for (t1,q1),(t2,q2)in zip(zip(co[0],co[1]),zip(co[0][1:],co[1][1:])):
    dt=t2-t1;dq=q2-q1
    if dt and dq:
     assert dt==dq
     for i in range(dt):
      paired.append((t1+i,q1+i))
      if ref['sequence'][t1+i]!=query[q1+i]:
       pos=q1+i;oldpos=pos if ori=='+'else len(q)-pos-1
       mm.append({'target_position0':t1+i,'oriented_query_position0':pos,'stored_query_position0':oldpos,'reference_base':ref['sequence'][t1+i],'query_base':query[pos],'BQ_by_member':[{'name_sha256':m['name_sha256'],'read1':m['read1'],'BQ':m['qualities'][oldpos]}for m in members]})
    elif dq:ins.append({'target_position0':t1,'query_start0':q1,'query_end0':q2,'sequence':query[q1:q2],'at_physical_reference_end':t1 in(0,len(ref['sequence']))})
    elif dt:dels.append({'target_start0':t1,'target_end0':t2,'query_position0':q1,'reference_end_overhang':q1 in(0,len(query))})
   r={'query_group':qi,'query_orientation':ori,'query_sequence':query,'query_length':len(query),'names':sorted({m['name_sha256']for m in members}),'primary_alignment_members':members,'reference':ref,'score':al.score,'coordinates':co,'alignment':str(al),'mismatches':mm,'query_insertions':ins,'reference_deletions':dels,'aligned_query_bases':len(paired),'aligned_query_fraction':len(paired)/len(query),'matching_bases':len(paired)-len(mm),'left_join_aligned_bases':sum(t<250 for t,q in paired)if ref['id'].startswith('proposed')else None,'right_join_aligned_bases':sum(t>=250 for t,q in paired)if ref['id'].startswith('proposed')else None}
   results.append(r);export.append({'case_id':f'DNA_Q{qi}{ori}_{ref["id"]}','query_sequence':query,'target_sequence':ref['sequence'],'reference':{k:v for k,v in ref.items()if k!='sequence'},'score':al.score,'coordinates':co})
report={'scoring':str(a),'bounded_scope':'Same frozen score conventions; two contiguous 750nt parent windows and the already validated proposed750nt model, both orientations; not a new genome-wide search.','source_records':rows,'distinct_primary_stored_sequences':len(group),'name_count':len(names),'results':results}
(P/'linked-read-alignments.json').write_text(json.dumps(report,indent=2)+'\n');(P/'linked-root-rescore-cases.json').write_text(json.dumps(export,indent=2)+'\n')
for r in results:print(r['query_group'],r['query_orientation'],r['reference']['id'],r['score'],r['coordinates'],r['alignment'])
