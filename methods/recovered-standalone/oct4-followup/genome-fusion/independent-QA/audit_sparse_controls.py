from pathlib import Path
from collections import defaultdict
import json,pysam
ROOT=Path(__file__).resolve().parents[1];P=ROOT/'compact-controls';all_records={}
for d in ('D1','D8'):
 out=defaultdict(list)
 with pysam.AlignmentFile(str(P/d/'Aligned.control.bam'),'rb') as f:
  for r in f:
   out[r.query_name].append({'flag':r.flag,'ref':r.reference_name,'start0':r.reference_start,'end0':r.reference_end,'mapq':r.mapping_quality,'cigar':r.cigarstring,'mate_ref':r.next_reference_name,'mate_start0':r.next_reference_start,'tlen':r.template_length,'tags':dict(r.tags),'query_length':r.query_length})
 all_records[d]=dict(out)
changed={n:{d:all_records[d].get(n) for d in all_records} for n in set(all_records['D1'])|set(all_records['D8']) if all_records['D1'].get(n)!=all_records['D8'].get(n)}
out={'changed_names':changed,'junctions':{d:[l.strip() for l in (P/d/'STAR.Chimeric.out.junction').read_text().splitlines() if 'BCR-ABL1-28\t' in l] for d in all_records}}
(Path(__file__).parent/'sparse-control-read-differences.json').write_text(json.dumps(out,indent=2))
