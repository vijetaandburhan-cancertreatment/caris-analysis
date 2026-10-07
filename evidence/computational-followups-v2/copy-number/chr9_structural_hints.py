"""Bounded discordant-pair/split-alignment hints, not structural-variant calls."""
import pathlib,json,collections,hashlib,re,csv,time
import pysam
OUT=pathlib.Path(__file__).resolve().parent;ROOT=OUT.parents[2];SRC=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
bam=pysam.AlignmentFile(str(SRC/'DNA_TN26-279853.bam'),'rb',index_filename=str(ROOT/'work/oct1-analysis/DNA_TN26-279853.bam.bai'))
windows=[('chr9',18000000,35000000),('chr3',51400000,53400000),('chr5',86300000,88300000)];pairs=collections.defaultdict(dict);splits=collections.defaultdict(dict);clips=collections.defaultdict(dict);qc=collections.Counter();started=time.time()
for ch,lo,hi in windows:
 for r in bam.fetch(ch,lo,hi):
  qc[ch+'_records_seen']+=1
  if r.flag&(4|256|512|1024|2048) or r.mapping_quality<30:continue
  rg=r.get_tag('RG') if r.has_tag('RG') else '';key=rg+'|'+r.query_name
  if r.is_read1 and not r.mate_is_unmapped and r.next_reference_id==r.reference_id and abs(r.next_reference_start-r.reference_start)>=1000:
   a,b=sorted([r.reference_start,r.next_reference_start]);orient=('R' if r.is_reverse else 'F')+('R' if r.mate_is_reverse else 'F')
   # Read1 orientation is not fixed left-to-right; use ordered orientation.
   if r.reference_start>r.next_reference_start:orient=orient[::-1]
   cluster=(ch,a//500,b//500,orient);pairs[cluster][key]={'left_start0':a,'right_start0':b,'read1_MAPQ':r.mapping_quality,'mate_MAPQ':r.get_tag('MQ') if r.has_tag('MQ') else None,'TLEN':r.template_length,'proper_pair':r.is_proper_pair,'cigar':r.cigarstring}
  if r.has_tag('SA'):
   for sa in r.get_tag('SA').rstrip(';').split(';'):
    z=sa.split(',')
    if len(z)<6 or int(z[4])<30:continue
    sch,sp,strand,scig,smq,nm=z[:6];sp=int(sp)-1
    if sch!=ch or abs(sp-r.reference_start)<1000:continue
    a,b=sorted([r.reference_start,sp]);cluster=(ch,a//500,b//500,('R' if r.is_reverse else 'F')+strand)
    splits[cluster][key]={'primary_start0':r.reference_start,'primary_end0':r.reference_end,'primary_cigar':r.cigarstring,'primary_MAPQ':r.mapping_quality,'supplementary_start0':sp,'supplementary_cigar':scig,'supplementary_MAPQ':int(smq),'supplementary_NM':int(nm)}
  for idx,(op,n) in enumerate(r.cigartuples or []):
   if op==4 and n>=20:
    endpoint=r.reference_start if idx==0 else r.reference_end;side='left' if idx==0 else 'right';cluster=(ch,endpoint//100,side);clips[cluster][key]={'endpoint0':endpoint,'clip_bases':n,'MAPQ':r.mapping_quality,'reverse':r.is_reverse,'primary_start0':r.reference_start,'cigar':r.cigarstring}
def pack(d,minn):
 rows=[]
 for cluster,vals in sorted(d.items(),key=lambda z:len(z[1]),reverse=True):
  if len(vals)<minn:continue
  rows.append({'cluster':list(cluster),'query_name_fragments':len(vals),'unique_alignment_descriptors':len({json.dumps(v,sort_keys=True) for v in vals.values()}),'fragment_evidence':[{'fragment_identifier_sha256':hashlib.sha256(k.encode()).hexdigest(),**v} for k,v in vals.items()]})
 return rows
result={'windows0based':windows,'elapsed_seconds':time.time()-started,'qc':dict(qc),'filters':'Primary mapped alignment MAPQ>=30, exclude secondary/supplementary/QCfail/duplicate. Samechrom mate/supplementary >=1kb away. Read1 only forpair clusters. SAalignmentMAPQ>=30. Groupstartcoordinates500bp;clipendpoint100bp.','pair_clusters_ge3':pack(pairs,3),'split_clusters_ge2':pack(splits,2),'clip_clusters_ge5':pack(clips,5),'limitations':['DiscordantpairmateMAPQunavailablewhenMQtagabsent. Longalignmentseparationisnotitselfproof ofdeletion.','Clustersbyalignmentstarts—not exactbreakpoint—canfragmentoneevent ormerge differentones.','SA/softclips require exactreferencealignment/assemblyandmappabilityaudit before structuralcall.','Lack ofclusters in exome/targetenricheddata doesnotexcludeanevent.','No somatic status or absolute copy inference.']}
(OUT/'regional-structural-hints.json').write_text(json.dumps(result,indent=2)+'\n')
print('QC',dict(qc),'pairs',len(result['pair_clusters_ge3']),'splits',len(result['split_clusters_ge2']),'clips',len(result['clip_clusters_ge5']))
for kind in ['pair_clusters_ge3','split_clusters_ge2']:
 print(kind)
 for r in result[kind][:20]:print(r['cluster'],r['query_name_fragments'],r['unique_alignment_descriptors'])
