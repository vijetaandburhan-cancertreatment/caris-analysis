#!/usr/bin/env python3
"""Separate selected-DNA sequence-pattern recount against independently supplied expected counts."""
import collections,hashlib,json,pathlib,pysam
base=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis');ws=pathlib.Path(__file__).resolve().parent;a=json.load(open(base/'oct4-fusion-read-audit/patient/junction-read-audit.json'));expect=json.load(open(ws/'independent-priority-review/bounded-DNA-details.json'));out=[]
bam=pysam.AlignmentFile(str(base/'TN26-279853/DNA_TN26-279853.bam'),'rb',index_filename=str(ws.parents[2]/'oct1-analysis/DNA_TN26-279853.bam.bai'))
# Explicit paths avoid depending on source caller implementation.
for cid in ['078bc453229c','e2a8e3625460']:
 r=next(x for x in a['rows'] if x['candidate_id']==cid);ch,p1=r['breakpoint1'].split(':');p1=int(p1);p2=int(r['breakpoint2'].split(':')[1]);left,right=r['fusion_transcript'].upper().split('|');intervals=[(p2-501,p1+500)] if p1-p2<1000 else [(p1-501,p1+500),(p2-501,p2+500)];patterns={arm:left[-arm:]+right[:arm] for arm in [25,50]};fragments=collections.defaultdict(set);alignments=set();n=0
 for start,end in intervals:
  for read in bam.fetch(ch,start,end):
   identity=(read.query_name,read.flag,read.reference_start,read.cigarstring)
   if identity in alignments:continue
   alignments.add(identity)
   if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail or read.is_duplicate or read.mapping_quality<20 or not read.query_sequence or read.query_qualities is None:continue
   n+=1;key=(read.query_name,read.get_tag('RG') if read.has_tag('RG') else '')
   for arm,p in patterns.items():
    hits=[];start=0
    while True:
     at=read.query_sequence.find(p,start)
     if at<0:break
     hits.append(min(read.query_qualities[at:at+len(p)]));start=at+1
    for q in [20,30]:
     if any(v>=q for v in hits):fragments[(arm,q)].add(key)
 z={'candidate_id':cid,'gene':r['gene1'],'eligible_read_alignments':n,'counts':[]}
 for arm in [25,50]:
  for q in [20,30]:
   v=len(fragments[(arm,q)]);e=expect[cid]['longer_exact_junction_markers'][f'{arm}+{arm}_Q{q}']['fragments'];z['counts'].append({'arm':arm,'minimum_base_quality':q,'separate_recount_fragments':v,'expected':e,'match':v==e})
 out.append(z)
r={'status':'passed' if all(z['match'] for row in out for z in row['counts']) else 'mismatch','rows':out,'method':'Fresh indexed DNA BAM fetch; explicit MAPQ20 and primary/mapped/nonduplicate/non-QCfail flag conditions; all string occurrences checked, Q20/Q30, query-name+read-group collapse. Same predeclared local intervals as peer; no CIGAR/haplotype or protein inference in this spot test.','limits':['Confirms local exact sequence evidence, not somatic status, genomic uniqueness, original molecule count or clinical relevance.','This is a second implementation over the same BAM and caller-derived marker, not independent tissue or assay validation.']};(ws/'independent-DNA-spot-recount.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
