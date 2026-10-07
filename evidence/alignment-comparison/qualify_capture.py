from pathlib import Path
import pysam,json,hashlib
from select_capture import capture
R=Path(__file__).resolve().parent;D=R/'capture-qualification';D.mkdir(exist_ok=True)
h={'HD':{'VN':'1.6','SO':'unsorted'},'SQ':[{'SN':'chr1','LN':1000},{'SN':'chr2','LN':1000}]}
rows=[]
def make(name,rid,start,cigar,flag=0,mq=0):
 r=pysam.AlignedSegment();r.query_name=name;r.reference_id=rid;r.reference_start=start;r.cigarstring=cigar;r.flag=flag;r.mapping_quality=mq
 n=sum(n for op,n in (r.cigartuples or []) if op in(0,1,4,7,8)) if cigar else 60
 r.query_sequence='A'*n;r.query_qualities=[35]*n;r.set_tag('ZZ','preserve');return r
rows=[make('site-only',0,70,'60M',mq=0),make('named',1,500,'60M',flag=2048,mq=0),make('named',-1,-1,None,flag=4),make('named',0,70,'60M',flag=256),make('skip-not-named',0,70,'10M100N10M'),make('outside',1,70,'60M'),make('deletion-only',0,70,'10M40D10M')]
with pysam.AlignmentFile(str(D/'input.bam'),'wb',header=h) as f:
 for r in rows:f.write(r)
(D/'sites.bed').write_text('chr1\t100\t101\n');(D/'names.txt').write_text('named\n')
capture(D/'input.bam',D/'sites.bed',D/'names.txt',D/'capture.bam',D/'result.json')
with pysam.AlignmentFile(str(D/'input.bam'),'rb') as f:a=[r.to_string() for r in f]
with pysam.AlignmentFile(str(D/'capture.bam'),'rb') as f:b=[r.to_string() for r in f]
assert b==a[:4],(len(a),len(b));x={'status':'PASS','input_records':7,'expected_retained':4,'exact_SAM_field_strings_preserved':True,'qualifies':['site overlap retained withMAPQ0','named mate atotherchromosome retained','namedunmapped/supplementary/secondary retained','recordmeetingbothruleswrittenonce','unnamedskip-only/deletion-only/outside omitted'],'scope':'Execution plumbing only; no clinical sensitivity claim'};(D/'qualification.json').write_text(json.dumps(x,indent=2)+'\n');print('Capture qualification PASS')
