from pathlib import Path
import json,collections,time,resource,re,heapq
import pysam
from Bio import Align
P=Path(__file__).parent;B=P.parent.parent.parent
N=json.loads((P/'nominations.json').read_text());Q={x['query']['oriented_id']:x for x in N['queries']}
fa=pysam.FastaFile(str(B/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'))
a=Align.PairwiseAligner(mode='global',match_score=2,mismatch_score=-4)
for typ in ('insertion','deletion'):
 setattr(a,f'open_internal_{typ}_score',-6);setattr(a,f'extend_internal_{typ}_score',-1)
for side in ('left','right'):
 setattr(a,f'open_{side}_deletion_score',0);setattr(a,f'extend_{side}_deletion_score',0);setattr(a,f'open_{side}_insertion_score',-2);setattr(a,f'extend_{side}_insertion_score',-2)
assert a.score('ACGT','TTACGT')==4 and a.score('TTACGTAA','ACGT')==8 and a.score('GGACGTGG','TTACGT')==1
best=collections.defaultdict(list);home=collections.defaultdict(list);ns=collections.Counter();hist=collections.defaultdict(collections.Counter);t0=time.time()
def offer(oid,r):
 q=Q[oid]['query'];r['score']=a.score(r['sequence'],q['oriented_sequence']);ns[oid]+=1;hist[oid][r['score']]+=1
 best[oid].append(r);best[oid].sort(key=lambda x:(-x['score'],x['kind'],x['id'],x['start0']));best[oid]=best[oid][:8]
 if q['kind']=='patient_representative':
  site=q['site'];is_home=(r['kind']=='genome' and r['id']==site['chrom'] and r['start0']<=site['pos1']-1<r['end0'])or(r['kind']=='transcript' and r.get('gene_name') in site['gene_symbols'].split(';'))
  if is_home:
   home[oid].append(r);home[oid].sort(key=lambda x:(-x['score'],x['kind'],x['id'],x['start0']));home[oid]=home[oid][:5]
for oid,v in Q.items():
 for w in v['references']['genome']:
  s=fa.fetch(w['reference_id'],w['start0'],w['end0']).upper()
  if s:offer(oid,{'kind':'genome','id':w['reference_id'],'start0':w['start0'],'end0':w['start0']+len(s),'sequence':s})
bytx=collections.defaultdict(list)
for oid,v in Q.items():
 for w in v['references']['transcript']:bytx[w['reference_id']].append((oid,w))
with pysam.FastxFile(str(B/'oct4-hla-allele-support/gencode.v37.transcripts.fa.gz'))as f:
 for tx in f:
  fields=tx.name.split('|');tid=fields[0]
  if tid not in bytx:continue
  for oid,w in bytx[tid]:
   s=tx.sequence[w['start0']:w['end0']].upper()
   offer(oid,{'kind':'transcript','id':tid,'gene_id':fields[1],'gene_name':fields[5] if len(fields)>5 else '?','transcript_name':fields[4] if len(fields)>4 else '?','start0':w['start0'],'end0':w['start0']+len(s),'sequence':s})
# A reference sequence following each original CIGAR splice path is an explicit comparator,
# not proof that this alignment, allele, or splice path is biologically correct.
paths=collections.defaultdict(list)
for oid,v in Q.items():
 q=v['query']
 if q['kind']!='patient_representative':continue
 seen=set()
 for m in q['members']:
  st=m['start0'];rp=st;parts=[fa.fetch(q['site']['chrom'],max(0,st-20),st).upper()];offset=len(parts[0]);siteindex=None
  for length,op in re.findall(r'(\d+)([MIDNSHP=X])',m['original_CIGAR']):
   n=int(length)
   if op in 'M=XD':
    if rp<=q['site']['pos1']-1<rp+n:siteindex=offset+q['site']['pos1']-1-rp
    s=fa.fetch(q['site']['chrom'],rp,rp+n).upper();parts.append(s);offset+=len(s);rp+=n
   elif op=='N':rp+=n
  parts.append(fa.fetch(q['site']['chrom'],rp,rp+20).upper());s=''.join(parts)
  for label,base in [('public_reference',None),('DNA_expected_allele',q['site']['expected_base']),('RNA_unexpected_allele',q['site']['unexpected_base'])]:
   ss=s if base is None or siteindex is None else s[:siteindex]+base+s[siteindex+1:]
   key=(ss,siteindex,label)
   if key in seen:continue
   seen.add(key);paths[oid].append({'kind':'source_CIGAR_reference_path','id':f'{q["site"]["chrom"]}:{st}:{m["original_CIGAR"]}:{label}','sequence':ss,'start0':0,'end0':len(ss),'site_index0':siteindex,'allele_model':label,'score':a.score(ss,q['oriented_sequence'])})
 paths[oid].sort(key=lambda r:-r['score'])

def detail(r,q):
 al=a.align(r['sequence'],q['oriented_sequence'])[0];co=al.coordinates.tolist();mm=[];ins=[];dels=[];paired=[]
 for (t1,q1),(t2,q2)in zip(zip(co[0],co[1]),zip(co[0][1:],co[1][1:])):
  dt=t2-t1;dq=q2-q1
  if dt and dq:
   assert dt==dq
   for i in range(dt):
    paired.append((t1+i,q1+i))
    if r['sequence'][t1+i]!=q['oriented_sequence'][q1+i]:mm.append({'target0':t1+i,'query0':q1+i,'reference':r['sequence'][t1+i],'query':q['oriented_sequence'][q1+i]})
  elif dq:ins.append([t1,q1,q2])
  elif dt:dels.append([t1,t2,q1,q1 in(0,len(q['oriented_sequence']))])
 return dict(r,alignment={'score':al.score,'coordinates':co,'text':str(al),'query_length':len(q['oriented_sequence']),'query_bases_aligned':len(paired),'aligned_fraction':len(paired)/len(q['oriented_sequence']),'mismatches':mm,'insertions':ins,'deletions':dels})
results=[];exports=[]
for oid,v in Q.items():
 q=v['query'];r={'query':q,'windows_scored':ns[oid],'best_score_window_count':hist[oid][best[oid][0]['score']]if best[oid]else 0,'top_alternatives':[detail(x,q)for x in best[oid]],'top_annotated_home':[detail(x,q)for x in home[oid]],'source_CIGAR_paths':[detail(x,q)for x in paths[oid][:3]],'source_CIGAR_path_count':len(paths[oid])}
 results.append(r)
 for label,ls in [('normal',r['top_alternatives'][:2]),('source_path',r['source_CIGAR_paths'][:1])]:
  for i,x in enumerate(ls):exports.append({'case_id':f'{oid}_{label}_{i}','query_id':q['query_id'],'query_sequence':q['oriented_sequence'],'target_sequence':x['sequence'],'score':x['score'],'reference':{k:v for k,v in x.items()if k not in('sequence','alignment')},'coordinates':x['alignment']['coordinates']})
(P/'scored-results.json').write_text(json.dumps({'results':results,'scored_windows':sum(ns.values()),'elapsed_seconds':time.time()-t0,'peak_RSS':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'aligner':str(a),'seed_caps':N['seed_caps'],'window_caps':N['window_caps'],'limitations':'Representative full reads only; caps/seed nomination not exhaustive; source-CIGAR reference paths and allele substitutions are comparators, not validated origins.'},indent=2)+'\n');(P/'independent-rescore-cases.json').write_text(json.dumps(exports,indent=2)+'\n')
print('scored',sum(ns.values()),'seconds',time.time()-t0,'RSS',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
