#!/usr/bin/env python3
import pathlib,json,pysam,collections,hashlib,resource,shutil,datetime
O=pathlib.Path(__file__).parent;P=O.parent;B=pathlib.Path.home()/'.local/share/codex/caris-analysis';DNA=B/'TN26-279853/DNA_TN26-279853.bam';INDEX=pathlib.Path.cwd()/'work/oct1-analysis/DNA_TN26-279853.bam.bai';fa=pysam.FastaFile(str(B/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'));models=json.loads((P/'proposed-models.json').read_text());assert models[0]['sequence']==models[1]['sequence'];model=models[0]['sequence'];join=models[0]['junction_index0'];rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1];sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
priorRNA=sha(P/'checksums.json');before=DNA.stat();assert shutil.disk_usage(O).free>=3*2**30
bps={'CCND1':{'chrom':'chr11','positions1':list(range(69651217,69651235)),'retained_side':'higher'},'HEG1':{'chrom':'chr3','positions1':list(range(125013572,125013590)),'retained_side':'lower'}}
for x in bps.values():x['fetch0']=[min(x['positions1'])-1-1000,max(x['positions1'])+1000]
markers=collections.defaultdict(list)
for delta in range(-17,1):
 j=join+delta
 for n in(50,100):
  seq=model[j-n//2:j+n//2];assert len(seq)==n
  for ori,s in [('+',seq),('-',rc(seq))]:markers[s].append({'length':n,'delta':delta,'orientation':ori})
assert len(markers)==72
# Independent synthetic pattern/quality logic regression, before private BAM.
for seq,roles in markers.items():
 assert seq in('AAA'+seq+'TTT') and min([30]*len(seq))>=30
 assert min([30]*(len(seq)-1)+[19])<20
bam=pysam.AlignmentFile(str(DNA),'rb',index_filename=str(INDEX));matebam=pysam.AlignmentFile(str(DNA),'rb',index_filename=str(INDEX));records={};regions_counts={}
def rec_key(a):return(a.query_name,a.flag,a.reference_id,a.reference_start,a.cigarstring,a.next_reference_id,a.next_reference_start,a.query_sequence)
def frag_key(a):return(a.query_name,a.get_tag('RG')if a.has_tag('RG')else'')
def strict(a):return not(a.is_unmapped or a.is_secondary or a.is_supplementary or a.is_duplicate or a.is_qcfail)and a.mapping_quality>=20
for gene,g in bps.items():
 nn=0
 for a in bam.fetch(g['chrom'],*g['fetch0']):records[rec_key(a)]=a;nn+=1
 regions_counts[gene]=nn
initial_record_count=len(records)

def exact(a):
 if not a.query_sequence or a.query_qualities is None:return[]
 out=[]
 for seq,roles in markers.items():
  i=a.query_sequence.find(seq)
  while i>=0:
   qq=min(a.query_qualities[i:i+len(seq)])
   if qq>=20:out.extend({'length':x['length'],'delta':x['delta'],'orientation':x['orientation'],'start0':i,'end0':i+len(seq),'BQmin':qq}for x in roles)
   i=a.query_sequence.find(seq,i+1)
 return out

def geometry(a):
 own=[(gene,g)for gene,g in bps.items()if not a.is_unmapped and a.reference_name==g['chrom']and a.reference_start<g['fetch0'][1]and a.reference_end>g['fetch0'][0]]
 info=[]
 for gene,g in own:
  other=next(v for k,v in bps.items()if k!=gene);record={'gene':gene,'expected_boundary_softclip':False,'other_softclip_near_boundary':False,'mate_other_locus':False,'SA_other_locus':False,'clips':[],'SA_links':[]}
  if a.is_paired and not a.mate_is_unmapped and a.next_reference_id>=0 and a.next_reference_name==other['chrom']and other['fetch0'][0]<=a.next_reference_start<other['fetch0'][1]:record['mate_other_locus']=True
  cig=a.cigartuples or[]
  for side,op in [('left',cig[0]if cig else(-1,0)),('right',cig[-1]if cig else(-1,0))]:
   if op[0]!=4 or op[1]<12:continue
   boundary=a.reference_start if side=='left'else a.reference_end
   expected=[p-1 if side=='left'else p for p in g['positions1']];distance=min(abs(boundary-p)for p in expected)
   if distance<=100:
    qs=0 if side=='left'else len(a.query_sequence)-op[1];clipped=a.query_sequence[qs:qs+op[1]];quality=a.query_qualities[qs:qs+op[1]]if a.query_qualities else[]
    record['clips'].append({'side':side,'length':op[1],'boundary0':boundary,'min_distance_to_equivalent_boundary':distance,'sequence':clipped,'BQmin':min(quality)if quality else None})
    if distance<=35 and((gene=='CCND1'and side=='left')or(gene=='HEG1'and side=='right')):record['expected_boundary_softclip']=True
    else:record['other_softclip_near_boundary']=True
  if a.has_tag('SA'):
   for ent in a.get_tag('SA').split(';'):
    if not ent:continue
    c,p,strand,cigar,mapq,nm=ent.split(',');pos=int(p)-1
    if c==other['chrom']and other['fetch0'][0]<=pos<other['fetch0'][1]:record['SA_other_locus']=True;record['SA_links'].append({'chrom':c,'start0':pos,'strand':strand,'CIGAR':cigar,'MAPQ':int(mapq),'NM':int(nm)})
  if any(record[k]for k in('expected_boundary_softclip','other_softclip_near_boundary','mate_other_locus','SA_other_locus')):info.append(record)
 return info
# Prioritize mate requests, collapse names+RG. All already-fetched mates remain in records.
candidates={}
for a in list(records.values()):
 hh=exact(a);gg=geometry(a);priority=0
 if hh:priority=4
 if any(x['SA_other_locus']or x['mate_other_locus']for x in gg):priority=max(priority,3)
 if any(x['expected_boundary_softclip']for x in gg):priority=max(priority,2)
 if gg:priority=max(priority,1)
 if priority and a.is_paired:
  k=frag_key(a)
  if k not in candidates or priority>candidates[k][0]:candidates[k]=(priority,a)
mate_lookups=[];selected=sorted(candidates.items(),key=lambda x:(-x[1][0],x[0]))[:500]
for key,(priority,a)in selected:
 if any(frag_key(z)==key and z.is_read1!=a.is_read1 for z in records.values()):mate_lookups.append({'name_sha256':hashlib.sha256(key[0].encode()).hexdigest(),'priority':priority,'status':'already_fetched'});continue
 try:
  mate=matebam.mate(a)
  if frag_key(mate)!=key:raise ValueError('mate name/RG does not match')
  records[rec_key(mate)]=mate;status='mate_added'
 except (ValueError,OSError)as e:status='unavailable:'+str(e)
 mate_lookups.append({'name_sha256':hashlib.sha256(key[0].encode()).hexdigest(),'priority':priority,'status':status})
# Recount all saved records, no duplicates by stored alignment descriptor.
matchsets=collections.defaultdict(set);hitrecords=[];flagcounts=collections.Counter();geomout=[];geocount=collections.defaultdict(set)
for a in records.values():
 flags=[]
 for yes,label in[(a.is_unmapped,'unmapped'),(a.is_secondary,'secondary'),(a.is_supplementary,'supplementary'),(a.is_duplicate,'duplicate'),(a.is_qcfail,'QCfail'),(a.mapping_quality<20,'MAPQ_lt20')]:
  if yes:flagcounts[label]+=1;flags.append(label)
 hh=exact(a)
 if hh:
  hitrecords.append({'name_sha256':hashlib.sha256(a.query_name.encode()).hexdigest(),'RG':frag_key(a)[1],'flag':a.flag,'MAPQ':a.mapping_quality,'reference':a.reference_name if not a.is_unmapped else None,'position1':a.reference_start+1,'CIGAR':a.cigarstring,'strict':strict(a),'flags':flags,'sequence':a.query_sequence,'qualities':list(a.query_qualities),'hits':hh})
  for n in(50,100):
   for q in(20,30):
    if any(h['length']==n and h['BQmin']>=q for h in hh):
     matchsets[f'all_flags_{n}nt_BQ{q}'].add(frag_key(a))
     if strict(a):matchsets[f'strict_{n}nt_BQ{q}'].add(frag_key(a))
     for flag in flags:matchsets[f'flag_{flag}_{n}nt_BQ{q}'].add(frag_key(a))
 gg=geometry(a)
 if gg:
  z={'name_sha256':hashlib.sha256(a.query_name.encode()).hexdigest(),'flag':a.flag,'MAPQ':a.mapping_quality,'chrom':a.reference_name,'pos1':a.reference_start+1,'CIGAR':a.cigarstring,'strict':strict(a),'mate_chr':a.next_reference_name if a.next_reference_id>=0 else None,'mate_pos1':a.next_reference_start+1,'template_length':a.template_length,'geometry':gg}
  geomout.append(z)
  for g in gg:
   for label in('expected_boundary_softclip','other_softclip_near_boundary','mate_other_locus','SA_other_locus'):
    if g[label]:geocount[(g['gene'],label,'strict'if strict(a)else'relaxed_only')].add(frag_key(a))
# Coverage on both sides, including reference-exact and non-reference aligned anchor counts.
coverage={}
for gene,g in bps.items():
 relevant=[a for a in records.values()if strict(a)and not a.is_unmapped and a.reference_name==g['chrom']and a.query_qualities]
 cov={};anchor={}
 for p1 in g['positions1']:
  names=set();basecounts=collections.defaultdict(set)
  for a in relevant:
   for qp,rp in a.get_aligned_pairs(matches_only=True):
    if rp==p1-1 and a.query_qualities[qp]>=20:names.add(frag_key(a));basecounts[a.query_sequence[qp]].add(frag_key(a))
  cov[str(p1)]={'fragments_aligned_BQ20':len(names),'bases':{k:len(v)for k,v in basecounts.items()}}
  for n in(25,50):
   start=p1-1 if g['retained_side']=='higher'else p1-n;end=start+n;want=fa.fetch(g['chrom'],start,end).upper();aligned=set();refexact=set()
   for a in relevant:
    mp={rp:qp for qp,rp in a.get_aligned_pairs(matches_only=True)};ix=[mp.get(pos)for pos in range(start,end)]
    if None in ix or any(y!=x+1 for x,y in zip(ix,ix[1:]))or min(a.query_qualities[k]for k in ix)<20:continue
    aligned.add(frag_key(a))
    if ''.join(a.query_sequence[k]for k in ix)==want:refexact.add(frag_key(a))
   anchor[f'{p1}_{n}nt']={'start0':start,'end0':end,'aligned_contiguous_BQ20_fragments':len(aligned),'reference_exact_BQ20_fragments':len(refexact)}
 coverage[gene]={'bases_by_equivalent_breakpoint':cov,'retained_anchors':anchor,'scope':'Covered aligned retained/exonic anchor bases only; does not assess all possible intronic genomic breakpoints/capture gaps.'}
# Save the complete bounded alignment set for independent recount, not an indexed whole-source copy.
with pysam.AlignmentFile(str(O/'bounded-DNA-records.bam'),'wb',header=bam.header)as out:
 for a in records.values():out.write(a)
after=DNA.stat();assert before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns;assert sha(P/'checksums.json')==priorRNA
summary={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_DNA':str(DNA),'source_DNA_size':before.st_size,'source_DNA_mtime_ns':before.st_mtime_ns,'source_index':str(INDEX),'source_index_sha256':sha(INDEX),'original_DNA_unchanged_size_mtime':True,'prior_RNA_manifest_sha256':priorRNA,'regions':bps,'fetched_alignment_counts_by_region':regions_counts,'initial_distinct_alignment_records':initial_record_count,'final_distinct_alignment_records_including_mates':len(records),'final_distinct_name_RG_fragments':len({frag_key(a)for a in records.values()}),'strict_alignment_records':sum(strict(a)for a in records.values()),'excluded_flag_counts_nonexclusive':dict(flagcounts),'exact_marker_counts':{f'{scope}_{n}nt_BQ{q}':len(matchsets[f'{scope}_{n}nt_BQ{q}'])for scope in('strict','all_flags')for n in(50,100)for q in(20,30)},'flag_specific_marker_counts':{k:len(v)for k,v in matchsets.items()if k.startswith('flag_')},'marker_hit_alignments':hitrecords,'geometry_counts':{'|'.join(k):len(v)for k,v in geocount.items()},'mate_candidate_name_count':len(candidates),'mate_query_cap':500,'mate_query_cap_saturated':len(candidates)>500,'mate_queries':mate_lookups,'coverage':coverage,'marker_definitions':[{ 'sequence':s,'representations':rr}for s,rr in markers.items()],'marker_equivalent_representations_counted_once_per_name_RG':True,'synthetic_marker_quality_controls_pass':True,'bounded_BAM_sha256':sha(O/'bounded-DNA-records.bam'),'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'limits':['Indexed locus fetch plus selected indexed mates, not full-DNA or unmapped-both-mates search.','All18 RNA-breakpoint representations describe one joined-sequence hypothesis; they are collapsed, not independent events.','A negative does not exclude a genomic breakpoint elsewhere in an intron, outside capture, or differing from the RNA junction.','Aligned/exonic anchor coverage is not evidence that all relevant intronic DNA was captured.','No validated structural-variant sensitivity, somatic status or clinical classification is asserted.']}
(O/'dna-audit.json').write_text(json.dumps(summary,indent=2)+'\n');(O/'softclip-SA-mate-details.json').write_text(json.dumps(geomout,indent=2)+'\n');print(json.dumps({k:summary[k]for k in('fetched_alignment_counts_by_region','initial_distinct_alignment_records','final_distinct_alignment_records_including_mates','strict_alignment_records','exact_marker_counts','geometry_counts','mate_candidate_name_count','mate_query_cap_saturated','peak_RSS_bytes')}))
