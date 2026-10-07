from pathlib import Path
import json,collections,re
P=Path(__file__).parent;B=P.parent.parent
A=json.loads((P/'ordinary-reference-search/scored-results.json').read_text());selected={};txids=set()
for r in A['results']:
 if r['query']['kind']!='patient_representative':continue
 qid=r['query']['query_id']
 if r['top_alternatives']:
  z=r['top_alternatives'][0]
  if qid not in selected or z['score']>selected[qid]['hit']['score']:selected[qid]={'query':r['query'],'hit':z}
for v in selected.values():
 if v['hit']['kind']=='transcript':txids.add(v['hit']['id'])
exons=collections.defaultdict(list)
with(B/'genome-fusion-reference/gencode.v37.primary_assembly.annotation.gtf').open()as f:
 for l in f:
  if l.startswith('#'):continue
  q=l.rstrip().split('\t',8)
  if q[2]!='exon':continue
  m=re.search(r'transcript_id "([^"]+)"',q[8])
  if m and m.group(1)in txids:exons[m.group(1)].append((q[0],int(q[3])-1,int(q[4]),q[6]))
for tid,xx in exons.items():exons[tid]=sorted(set(xx),key=lambda x:x[1],reverse=xx[0][3]=='-')
def genomepos(tid,tp):
 for chrom,s,e,strand in exons[tid]:
  if tp<e-s:return chrom,s+tp if strand=='+'else e-1-tp,strand
  tp-=e-s
 return None
rows=[]
for qid,v in selected.items():
 q=v['query'];h=v['hit'];co=h['alignment']['coordinates'];source_positions=sorted({m['observation_qpos0']for m in q['members']});mapped=[]
 for sourceq in source_positions:
  qp=sourceq if q['orientation']=='+'else len(q['sequence'])-1-sourceq;ts=None;kind='query_insertion_or_terminal_gap'
  for (t1,q1),(t2,q2)in zip(zip(co[0],co[1]),zip(co[0][1:],co[1][1:])):
   if t2>t1 and q2>q1 and q1<=qp<q2:ts=t1+qp-q1;kind='aligned';break
  gp=None
  if ts is not None:
   if h['kind']=='genome':gp=(h['id'],h['start0']+ts,'+')
   else:gp=genomepos(h['id'],h['start0']+ts)
   kind='same_genomic_position'if gp and gp[0]==q['site']['chrom']and gp[1]==q['site']['pos1']-1 else 'different_genomic_position'
  mapped.append({'stored_BAM_query_position0':sourceq,'oriented_query_position0':qp,'classification':kind,'target_window_position0':ts,'genomic_mapping':gp,'target_base':h['sequence'][ts]if ts is not None else None,'query_base':q['oriented_sequence'][qp]})
 rows.append({'query_id':qid,'site':q['site'],'assay':q['assay'],'role':q['role'],'base':q['base'],'family_names':q['family_names'],'best_reference':{k:v for k,v in h.items()if k not in('sequence','alignment')},'orientation':q['orientation'],'mapping':mapped,'note':'One first-best affine alignment; repeat-equivalent traceback alternatives not enumerated, so placement ambiguity not a new genotype call.'})
(P/'candidate-position-mapping.json').write_text(json.dumps({'rows':rows,'transcript_exons':dict(exons),'limits':'The queried base may move or enter an insertion in a full-read optimal alignment. This does not establish biological origin or correct breakpoint/indel genotype.'},indent=2)+'\n')
for r in rows:
 if r['role']=='unexpected':print(r['query_id'],r['best_reference']['id'],r['best_reference']['score'],r['mapping'])
