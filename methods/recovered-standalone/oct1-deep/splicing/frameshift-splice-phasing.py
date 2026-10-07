"""Local physical linkage of reported indels to nearby annotated splice junctions."""
import collections,json,time
from pathlib import Path
import pysam
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
A=json.loads((OUT/'target-annotation.json').read_text())
BAM='/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853/RNA_TN26-279853.bam'; INDEX=str(ROOT/'work/oct1-analysis/RNA_TN26-279853.bam.bai')
results=[]
for gene,ds,dn in [('BAP1',52408550,11),('RASA1',87332560,1)]:
 g=A['genes'][gene];tx=A['transcripts'][g['selected_transcript']];ex=tx['exons'];canonical={(e,s) for (_,e),(s,_) in zip(ex,ex[1:])}
 frags=collections.defaultdict(lambda:{'alt':False,'junctions':set(),'direct_alt_junctions':set(),'examples':[]})
 with pysam.AlignmentFile(BAM,'rb',index_filename=INDEX) as b:
  for r in b.fetch(g['chrom'],g['start0'],g['end0']):
   if r.flag&(4|256|512|1024|2048) or r.mapping_quality<20 or r.has_tag('NH') and r.get_tag('NH')!=1:continue
   key=(r.get_tag('RG') if r.has_tag('RG') else '')+'|'+r.query_name
   pos=r.reference_start;qpos=0;alt=False;junctions=[];c=r.cigartuples;qual=r.query_qualities
   for i,(op,n) in enumerate(c):
    if op==2 and pos==ds and n==dn and 0<i<len(c)-1 and c[i-1][0] in (0,7,8) and c[i+1][0] in (0,7,8) and min(c[i-1][1],c[i+1][1])>=5 and qual is not None and min(qual[qpos-5:qpos+5])>=20:alt=True
    if op==3 and 0<i<len(c)-1 and c[i-1][0] in (0,7,8) and c[i+1][0] in (0,7,8) and min(c[i-1][1],c[i+1][1])>=12 and qual is not None and min(qual[qpos-12:qpos+12])>=20:junctions.append((pos,pos+n))
    if op in (0,2,3,7,8):pos+=n
    if op in (0,1,4,7,8):qpos+=n
   if not alt and not junctions:continue
   d=frags[key]; d['alt']|=alt;d['junctions'].update(junctions)
   if alt:d['direct_alt_junctions'].update(junctions)
   if alt:d['examples'].append({'query_name':r.query_name,'flag':r.flag,'start0':r.reference_start,'cigar':r.cigarstring,'junctions':junctions})
 by_junction=collections.defaultdict(lambda:{'paired_fragments':0,'direct_read_fragments':0,'examples':[]})
 for key,d in frags.items():
  if not d['alt']:continue
  for j in d['junctions']:
   x=by_junction[j];x['paired_fragments']+=1;x['direct_read_fragments']+=j in d['direct_alt_junctions']
   if len(x['examples'])<4:x['examples'].extend(d['examples'][:1])
 rows=[{'intron_start0':s,'intron_end0':e,'canonical_selected_transcript':(s,e) in canonical,**d} for (s,e),d in sorted(by_junction.items())]
 results.append({'gene':gene,'chrom':g['chrom'],'transcript':tx['transcript_id'],'genomic_deletion_start0':ds,'deleted_bases':dn,'qualifying_alt_fragment_names':sum(d['alt'] for d in frags.values()),'alt_fragment_names_with_any_qualified_junction':sum(d['alt'] and bool(d['junctions']) for d in frags.values()),'linked_junctions':rows})
report={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'method':'Primary MAPQ>=20,NH1; exact reported CIGAR deletion and >=5 matched BQ20 bases on each deletion side; junction>=12 immediately adjacent matched BQ20 bases each side. Pairing by RG+queryname. Direct-read counts require variant and junction in same alignment; paired-fragment counts may place them on mates. No UMI/full-length transcript/protein confirmation. Canonical means selected GENCODE37 MANE transcript junction membership, not full-isoform proof. These stricter conditions are not identical to prior indel evidence table filters.','results':results}
(OUT/'frameshift-splice-phasing.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({**report,'results':[{k:v for k,v in r.items() if k!='linked_junctions'}|{'linked_junctions':[{k:v for k,v in j.items() if k!='examples'} for j in r['linked_junctions']]} for r in results]},indent=2))
