#!/usr/bin/env python3
import pathlib,json,pysam,collections,time,resource,shutil,heapq,hashlib
from Bio import Align
O=pathlib.Path(__file__).parent;B=pathlib.Path.home()/'.local/share/codex/caris-analysis';N=json.loads((O/'nominations.json').read_text());models=json.loads((O/'proposed-models.json').read_text());fa=pysam.FastaFile(str(B/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'));t0=time.time();a=Align.PairwiseAligner(mode='global',match_score=2,mismatch_score=-4)
for gap in('insertion','deletion'):
 setattr(a,f'open_internal_{gap}_score',-6);setattr(a,f'extend_internal_{gap}_score',-1)
for side in('left','right'):
 setattr(a,f'open_{side}_deletion_score',0);setattr(a,f'extend_{side}_deletion_score',0);setattr(a,f'open_{side}_insertion_score',-2);setattr(a,f'extend_{side}_insertion_score',-2)
(O/'scoring-settings.txt').write_text(str(a)+'\nGap oflengthk=open+(k-1)*extension. Target=reference window; query=full oriented read.\n')
def guard():
 if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*2**30:raise RuntimeError('RSS4GiB')
 if shutil.disk_usage(O).free<3*2**30:raise RuntimeError('Free space3GiB')
# Independent hand-computable controls, not clinical sensitivity.
assert a.score('TTTACGTACCC','ACGTA')==10
assert a.score('ACGTAC','ACATAC')==6
assert a.score('ACGTACGTACTAG','ACGTAACGTACTAG')==20
queries={x['query']['oriented_id']:x for x in N['queries']};best={k:[]for k in queries};allcounts=collections.Counter();scorehist={k:collections.Counter()for k in queries};ctr=0

def offer(oid,refkind,refid,start,seq,metadata=None):
 global ctr
 q=queries[oid]['query']['oriented_sequence'];score=a.score(seq,q);ctr+=1;allcounts[oid]+=1;scorehist[oid][score]+=1
 item={'reference_kind':refkind,'reference_id':refid,'start0':start,'end0':start+len(seq),'sequence':seq,'score':score}
 if metadata:item.update(metadata)
 best[oid].append(item);best[oid].sort(key=lambda x:(-x['score'],x['reference_kind'],x['reference_id'],x['start0']));best[oid]=best[oid][:10]
 if ctr%20000==0:guard();print(json.dumps({'scored':ctr,'elapsed':time.time()-t0,'RSS':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}),flush=True)
for oid,v in queries.items():
 for w in v['references']['genome']:
  seq=fa.fetch(w['reference_id'],w['start0'],w['end0']).upper()
  if seq:offer(oid,'genome',w['reference_id'],w['start0'],seq)
print('genome scoring complete',ctr,flush=True)
bytx=collections.defaultdict(list)
for oid,v in queries.items():
 for w in v['references']['transcript']:bytx[w['reference_id']].append((oid,w))
found=set()
with pysam.FastxFile(str(B/'oct4-hla-allele-support/gencode.v37.transcripts.fa.gz'))as f:
 for tx in f:
  tid=tx.name.split('|')[0]
  if tid not in bytx:continue
  found.add(tid);s=tx.sequence.upper()
  for oid,w in bytx[tid]:
   seq=s[w['start0']:w['end0']]
   if seq:offer(oid,'transcript',tid,w['start0'],seq,{'transcript_name':tx.name})
assert found==set(bytx)
# Same semi-global scoring for proposed models, without pretending a parent-only alignment supports crossing.
proposed={}
for oid,v in queries.items():
 q=v['query']['oriented_sequence'];proposed[oid]=[{**m,'reference_kind':'proposed','reference_id':m['id'],'start0':0,'end0':len(m['sequence']),'score':a.score(m['sequence'],q)}for m in models]

def detailed(target,query,v):
 als=a.align(target,query);al=als[0];co=al.coordinates.tolist();mismatch=[];ins=[];deletion=[];pairs=[]
 for (t1,q1),(t2,q2)in zip(zip(co[0],co[1]),zip(co[0][1:],co[1][1:])):
  dt=t2-t1;dq=q2-q1
  if dt and dq:
   assert dt==dq
   for i in range(dt):
    pairs.append((t1+i,q1+i))
    if target[t1+i]!=query[q1+i]:mismatch.append({'query_position0':q1+i,'target_position0':t1+i,'reference':target[t1+i],'query':query[q1+i]})
  elif dq:ins.append({'query_start0':q1,'query_end0':q2,'target_position0':t1,'at_target_end':t1 in(0,len(target)),'bases':query[q1:q2]})
  elif dt:deletion.append({'target_start0':t1,'target_end0':t2,'query_position0':q1,'reference_terminal':q1 in(0,len(query))})
 mappedq=[qq for tt,qq in pairs];left=min(mappedq)if mappedq else len(query);right=len(query)-max(mappedq)-1 if mappedq else 0
 d={'score':al.score,'coordinates':co,'alignment_text':str(al),'query_length':len(query),'query_bases_aligned_to_reference':len(pairs),'aligned_query_fraction':len(pairs)/len(query),'query_terminal_unaligned_left':left,'query_terminal_unaligned_right':right,'query_internal_insertions':ins,'target_deletions_or_overhangs':deletion,'mismatches':mismatch,'matching_bases':len(pairs)-len(mismatch),'minimum_error_bases':len(mismatch)+sum(x['query_end0']-x['query_start0']for x in ins)+sum(x['target_end0']-x['target_start0']for x in deletion if not x['reference_terminal'])}
 # Qualities retain original mate identity and raw-read orientation.
 members=v['query'].get('members',[])
 for mm in mismatch:
  qo=mm['query_position0'];mm['base_qualities']=[{'pair_id':m['pair_id'],'mate':m['mate'],'raw_position0':qo if v['query']['orientation']=='+'else len(query)-1-qo,'quality':ord(m['quality'][qo if v['query']['orientation']=='+'else len(query)-1-qo])-33}for m in members]
 return d
out=[];export=[]
for oid,v in queries.items():
 query=v['query']['oriented_sequence'];vv={'query':v['query'],'nominated_windows_scored':allcounts[oid],'best_alternative_score':best[oid][0]['score']if best[oid]else None,'top_alternative_window_count_at_best_score':scorehist[oid][best[oid][0]['score']]if best[oid]else 0,'top_alternatives':[],'proposed_models':[]}
 for rank,w in enumerate(best[oid]):
  dd=detailed(w['sequence'],query,v);vv['top_alternatives'].append(dict(w,alignment=dd))
  if rank<3:export.append({'case_id':oid+'_alternative'+str(rank+1),'query_id':oid,'query_sequence':query,'target_sequence':w['sequence'],'reference':{k:w[k]for k in('reference_kind','reference_id','start0','end0')},'score':dd['score'],'coordinates':dd['coordinates']})
 for w in proposed[oid]:
  dd=detailed(w['sequence'],query,v);j=w['junction_index0'];co=dd['coordinates'];mapped=[]
  for (ts,qs),(te,qe)in zip(zip(co[0],co[1]),zip(co[0][1:],co[1][1:])):
   if te>ts and qe>qs:mapped.extend(range(ts,te))
  dd['mapped_bases_left_of_join']=sum(x<j for x in mapped);dd['mapped_bases_right_of_join']=sum(x>=j for x in mapped);dd['crosses_join_with12bp_each_side']=dd['mapped_bases_left_of_join']>=12 and dd['mapped_bases_right_of_join']>=12
  vv['proposed_models'].append(dict(w,alignment=dd));export.append({'case_id':oid+'_'+w['id'],'query_id':oid,'query_sequence':query,'target_sequence':w['sequence'],'reference':{k:w[k]for k in('reference_kind','reference_id','start0','end0')},'score':dd['score'],'coordinates':dd['coordinates']})
 out.append(vv)
report={'status':'complete','method':'Exact-seed nominations followed by equal semi-global affine scores; top windows are not certified unique global alignments','scored_windows':ctr,'elapsed_seconds':time.time()-t0,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'seed_caps':N['seed_caps'],'window_caps':N['window_caps'],'results':out}
(O/'scored-results.json').write_text(json.dumps(report,indent=2)+'\n');(O/'root-independent-rescore-cases.json').write_text(json.dumps(export,indent=2)+'\n')
with(O/'top-windows.fasta').open('w')as f:
 for x in export:f.write('>'+x['case_id']+' target\n'+x['target_sequence']+'\n>'+x['case_id']+' query\n'+x['query_sequence']+'\n')
print(json.dumps({k:report[k]for k in('status','scored_windows','elapsed_seconds','peak_RSS_bytes','seed_caps','window_caps')}))
