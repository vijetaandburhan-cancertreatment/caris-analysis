from pathlib import Path
import json,collections,gzip,hashlib,datetime,resource,shutil,csv,re
import pysam
P=Path(__file__).parent;B=P.parent.parent;WS=Path('/Users/burhanazeem/Documents/Codex/2026-09-05/finances-plugin-finances-openai-curated-remote-3')
S=json.loads((P.parent/'independent-count-audit/sensitivity-and-exact-discordances.json').read_text())
locs={}
for e in S['baseline_events_all_sensitivities']:
 locs[(e['chrom'],e['pos1'])]={'chrom':e['chrom'],'pos1':e['pos1'],'gene_symbols':e['gene_symbols'],'unexpected_base':e['unexpected_base'],'expected_base':max(e['DNA'],key=e['DNA'].get),'selection':'primary'}
for e in S['exact_events_by_mode']['MAPQ60']:
 k=(e['chrom'],e['pos1'])
 if k not in locs:locs[k]={'chrom':e['chrom'],'pos1':e['pos1'],'gene_symbols':e['gene_symbols'],'unexpected_base':e['unexpected_base'],'expected_base':max(e['DNA'],key=e['DNA'].get),'selection':'MAPQ60_only'}
assert len(locs)==9
fa=pysam.FastaFile(str(B/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'))
def obs(a,pos0):
 rp=a.reference_start;qp=0;hit=None;near=[]
 for op,n in a.cigartuples or []:
  if op in(0,7,8):
   if rp<=pos0<rp+n:hit={'qpos0':qp+pos0-rp,'matchblock_left':pos0-rp,'matchblock_right':rp+n-1-pos0,'block_start0':rp,'block_end0':rp+n}
   rp+=n;qp+=n
  elif op in(1,4):
   near.append({'op':'I'if op==1 else'S','length':n,'rpos0':rp,'qpos0':qp,'distance_genomic':abs(rp-pos0)});qp+=n
  elif op in(2,3):
   near.append({'op':'D'if op==2 else'N','length':n,'rpos0':rp,'rend0':rp+n,'qpos0':qp,'distance_genomic':min(abs(rp-pos0),abs(rp+n-pos0))});rp+=n
  elif op not in(5,6):raise RuntimeError(op)
 if hit:hit['near_CIGAR_events']=near
 return hit
# Coordinate/CIGAR tests, independent from primary pileup and prior direct-count code.
def rd(cigar,seq,start=100):
 x=pysam.AlignedSegment();x.reference_start=start;x.cigarstring=cigar;x.query_sequence=seq;x.query_qualities=[30]*len(seq);return x
assert obs(rd('20M','A'*20),100)['qpos0']==0
assert obs(rd('5S20M','C'*5+'A'*20),105)['qpos0']==10
assert obs(rd('10M100N10M','A'*20),210)['qpos0']==10
assert obs(rd('10M1D10M','A'*20),110)is None
assert obs(rd('10M1I10M','A'*21),110)['qpos0']==11
assert obs(rd('20M','A'*20),120)is None
(P/'synthetic-controls.json').write_text(json.dumps({'manual_CIGAR_position_tests':6,'pass':True,'note':'BQ/filter/fragment conflicts are also checked by independently reproducing baseline and clean counts; these CIGAR tests do not establish clinical sensitivity.'},indent=2)+'\n')
rows=[];summaries=[];reps=[]
for (chrom,pos1),site in sorted(locs.items(),key=lambda z:(int(z[0][0][3:]),z[0][1])):
 site=dict(site,reference_base=fa.fetch(chrom,pos1-1,pos1).upper(),genomic_context51=fa.fetch(chrom,pos1-26,pos1+25).upper())
 for assay in('DNA','RNA'):
  if shutil.disk_usage(P).free<3*2**30:raise RuntimeError('disk guard')
  per=[]
  with pysam.AlignmentFile(str(B/f'TN26-279853/{assay}_TN26-279853.bam'),'rb',index_filename=str(WS/f'work/oct1-analysis/{assay}_TN26-279853.bam.bai'))as bam:
   for a in bam.fetch(chrom,pos1-1,pos1):
    h=obs(a,pos1-1);seq=a.query_sequence;qual=a.query_qualities
    r={'site':f'{chrom}:{pos1}','assay':assay,'name':a.query_name,'RG':a.get_tag('RG')if a.has_tag('RG')else'','flag':a.flag,'MAPQ':a.mapping_quality,'NH':a.get_tag('NH')if a.has_tag('NH')else None,'NM':a.get_tag('NM')if a.has_tag('NM')else None,'SA':a.get_tag('SA')if a.has_tag('SA')else None,'chrom':a.reference_name,'start0':a.reference_start,'end0':a.reference_end,'CIGAR':a.cigarstring,'reverse':a.is_reverse,'read1':a.is_read1,'mate_chrom':a.next_reference_name,'mate_start0':a.next_reference_start,'TLEN':a.template_length,'sequence':seq,'qualities':list(qual)if qual is not None else None,'observation':h,'baseline':False,'clean':False}
    if h and seq and qual is not None:
     q=h['qpos0'];h.update(base=seq[q].upper(),BQ=qual[q],aligned_query_left=q-a.query_alignment_start,aligned_query_right=a.query_alignment_end-1-q,query_length=len(seq),query_context51=seq[max(0,q-25):q+26],context51_complete=q>=25 and len(seq)-q-1>=25)
     r['baseline']=not a.flag&0xF0C and bool(a.flag&2)and a.mapping_quality>=30 and qual[q]>=25 and h['base']in'ACGT'
     r['clean']=r['baseline']and not any(op==4 for op,n in a.cigartuples)and min(h['aligned_query_left'],h['aligned_query_right'])>=5
    per.append(r);rows.append(r)
  ss={'site':site,'assay':assay,'fetched_alignments':len(per),'modes':{},'allele_detail':{}}
  for mode in('baseline','clean'):
   names=collections.defaultdict(set)
   for r in per:
    if r[mode]:names[(r['name'],r['RG'])].add(r['observation']['base'])
   ct=collections.Counter(next(iter(v))for v in names.values()if len(v)==1);ss['modes'][mode]={'counts':{b:ct[b]for b in'ACGT'},'conflicting_names':sum(len(v)>1 for v in names.values())}
   if mode=='baseline':baseline_names=names
  for base in 'ACGT':
   rr=[r for r in per if r['baseline']and r['observation']['base']==base and baseline_names[(r['name'],r['RG'])]=={base}]
   if not rr:continue
   nameset={(r['name'],r['RG'])for r in rr};byseq=collections.defaultdict(list)
   for r in rr:byseq[r['sequence']].append(r)
   fams=[]
   for seq,xx in byseq.items():
    fams.append({'sequence':seq,'names':len({(r['name'],r['RG'])for r in xx}),'alignments':len(xx),'records':xx})
   fams.sort(key=lambda x:(-x['names'],x['sequence']))
   def hist(f):return dict(collections.Counter(str(f(r))for r in rr))
   ss['allele_detail'][base]={'names':len(nameset),'alignments':len(rr),'full_stored_sequence_families':len(fams),'distinct_alignment_starts':len({(r['start0'],r['CIGAR'],r['reverse'])for r in rr}),'hist':{'CIGAR':hist(lambda r:r['CIGAR']),'MAPQ':hist(lambda r:r['MAPQ']),'NH':hist(lambda r:r['NH']),'BQ':hist(lambda r:r['observation']['BQ']),'qpos0':hist(lambda r:r['observation']['qpos0']),'strand':hist(lambda r:r['reverse']),'has_softclip':hist(lambda r:any(e['op']=='S'for e in r['observation']['near_CIGAR_events'])),'matchblock_left':hist(lambda r:r['observation']['matchblock_left']),'matchblock_right':hist(lambda r:r['observation']['matchblock_right'])},'top_sequence_families':[{'sequence':x['sequence'],'names':x['names'],'alignments':x['alignments']}for x in fams[:12]]}
   limit=(8 if base==site['unexpected_base']else 4 if base==site['expected_base']else 0)if assay=='RNA'else 2 if base==site['unexpected_base']else 0
   for fi,f in enumerate(fams[:limit],1):
    reps.append({'id':f'{chrom}_{pos1}_{assay}_{base}_{fi}','site':site,'assay':assay,'base':base,'role':'unexpected'if base==site['unexpected_base']else'expected','sequence':f['sequence'],'family_names':f['names'],'represented_records':f['records']})
  # Expected exact baseline counts from independent whole-panel CIGAR audit, loaded separately.
  table=list(csv.DictReader((P.parent/f'independent-count-audit/{assay}-direct-counts.tsv').open(),delimiter='\t'));tr=next(r for r in table if r['chrom']==chrom and int(r['pos1'])==pos1)
  for mode,m in [('baseline','baseline'),('clean','clean_no_softclip_edge5')]:assert ss['modes'][mode]['counts']=={b:int(tr[f'{m}_{b}'])for b in'ACGT'},(site,assay,mode)
  summaries.append(ss);print(site['chrom'],site['pos1'],assay,ss['fetched_alignments'],ss['modes'],flush=True)
  assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<4*2**30
with gzip.open(P/'all-site-read-records.jsonl.gz','wt')as f:
 for r in rows:f.write(json.dumps(r,separators=(',',':'))+'\n')
(P/'read-context-summary.json').write_text(json.dumps({'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'loci':len(locs),'baseline_and_clean_counts_independently_match':True,'summaries':summaries,'records':len(rows),'peak_RSS':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},indent=2)+'\n');(P/'representative-reads.json').write_text(json.dumps(reps,indent=2)+'\n')
print('done',len(rows),'records',len(reps),'representatives')
