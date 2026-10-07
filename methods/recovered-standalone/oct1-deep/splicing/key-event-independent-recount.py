"""Independent pysam.find_introns check of exact specified exon-skipping junctions."""
from pathlib import Path
import json,time,collections
import pysam
OUT=Path(__file__).resolve().parent; ROOT=OUT.parents[2]
A=json.loads((OUT/'target-annotation.json').read_text()); GENES=A['genes']; TX=A['transcripts']
BAM='/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853/RNA_TN26-279853.bam'
INDEX=str(ROOT/'work/oct1-analysis/RNA_TN26-279853.bam.bai')
results=[]
for gene,exon_number in [('MET',14),('BAP1',4),('RASA1',3)]:
 g=GENES[gene];tx=TX[g['selected_transcript']];ex=tx['exons'];i=exon_number-1 if g['strand']=='+' else len(ex)-exon_number
 assert 0<i<len(ex)-1
 events={'skip':(ex[i-1][1],ex[i+1][0]),'left_inclusion':(ex[i-1][1],ex[i][0]),'right_inclusion':(ex[i][1],ex[i+1][0])}
 filters={'primary_flag_pass':0,'mapq_unique_pass':0}
 with pysam.AlignmentFile(BAM,'rb',index_filename=INDEX) as b:
  def iterator():
   for r in b.fetch(g['chrom'],ex[i-1][0],ex[i+1][1]):
    if r.flag & (4|256|512|1024|2048):continue
    filters['primary_flag_pass']+=1
    if r.mapping_quality<20 or r.has_tag('NH') and r.get_tag('NH')!=1:continue
    filters['mapq_unique_pass']+=1
    yield r
  counts=b.find_introns(iterator())
 r={'gene':gene,'transcript':tx['transcript_id'],'strand':g['strand'],'exon_number_transcript_1based':exon_number,'exon_start0':ex[i][0],'exon_end0':ex[i][1], 'filters':filters,'events':{name:{'chrom':g['chrom'],'intron_start0':s,'intron_end0':e,'read_count_no_anchor_filter':counts.get((s,e),0)} for name,(s,e) in events.items()}}
 results.append(r)
report={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'method':'pysam AlignmentFile.find_introns on independently indexed local fetch, primary flags excluding4/256/512/1024/2048; MAPQ>=20, NH=1 where present. No anchor-length/base-quality filters. Counts are alignments, not molecules. Exact zeros are broader-filter negatives for specified coordinates only, not absence of other splicing or low-frequency disease.','results':results}
(OUT/'key-event-independent-recount.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
