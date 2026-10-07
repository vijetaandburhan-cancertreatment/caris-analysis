"""Bounded DNA support check near two low-level junctions absent from queried GTEx atlas."""
from pathlib import Path
import csv,json,collections,time
import pysam
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[2]
rows=[r for r in csv.DictReader((OUT/'candidate-gtex-comparison.tsv').open(),delimiter='\t') if r['reference_aware_research_pass']=='1' and not r['GTEx_samples_with_junction']]
BAM='/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853/DNA_TN26-279853.bam';INDEX=str(ROOT/'work/oct1-analysis/DNA_TN26-279853.bam.bai')
prov=json.loads((OUT/'reference-motif-provenance.json').read_text());res=[]
for row in rows:
 gene=row['gene'];d=json.load(open(prov['references'][gene]['path']));dna=d['dna'].upper();base=d['start'];s=int(row['intron_start0']);e=int(row['intron_end0']);canonical=[int(x) for x in row['closest_selected_junction'].split('-')]
 positions=set()
 for x in [s,e,*canonical]:positions.update(range(x-12,x+12))
 seen=collections.defaultdict(lambda:collections.defaultdict(set)); strand=collections.defaultdict(lambda:collections.Counter());indels=collections.defaultdict(set)
 with pysam.AlignmentFile(BAM,'rb',index_filename=INDEX) as b:
  for r in b.fetch(row['chrom'],min(positions),max(positions)+1):
   if r.flag&(4|256|512|1024|2048) or r.mapping_quality<20:continue
   name=(r.get_tag('RG') if r.has_tag('RG') else '')+'|'+r.query_name
   pos=r.reference_start;q=0
   for op,n in r.cigartuples:
    if op in (0,7,8):
     for p in range(max(pos,min(positions)),min(pos+n,max(positions)+1)):
      qi=q+(p-pos)
      if p in positions and r.query_qualities[qi]>=20:
       letter=r.query_sequence[qi];seen[p][letter].add(name);strand[p][(letter,r.is_reverse)]+=1
    elif op in (1,2) and pos in positions:
     if q>=5 and q+5<=len(r.query_sequence) and min(r.query_qualities[q-5:q+5])>=20:
      indels[(pos,'I' if op==1 else 'D',n,r.query_sequence[q:q+n] if op==1 else dna[pos-base:pos-base+n])].add(name)
    if op in (0,2,3,7,8):pos+=n
    if op in (0,1,4,7,8):q+=n
 base_rows=[]
 for p in sorted(positions):
  reference=dna[p-base];bases=seen[p];base_rows.append({'position_1based':p+1,'reference':reference,'allele_fragment_counts':{a:len(v) for a,v in bases.items()},'allele_read_strands':{a:{'forward':strand[p][(a,False)],'reverse':strand[p][(a,True)]} for a in bases}})
 res.append({'gene':gene,'chrom':row['chrom'],'novel_intron_start0':s,'novel_intron_end0':e,'canonical_intron':canonical,'base_counts':base_rows,'indels':[{'start0':p,'type':typ,'length':n,'sequence':seq,'fragment_names':len(names)} for (p,typ,n,seq),names in sorted(indels.items())]})
report={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'method':'Existing DNA BAM primary MAPQ>=20 BQ>=20; flags unmapped/secondary/supplementary/QCfail/duplicate excluded. Targeted +/-12 bp around candidate and closest-canonical intron boundaries. Alleles collapsed by RG+name per base, without full haplotype or calibrated somatic calling. Indels require 5 BQ20 query bases on each side; adjacent matched CIGAR length not additionally required in this exploratory inventory. No sample-normal comparison.','results':res}
(OUT/'rare-junction-DNA-check.json').write_text(json.dumps(report,indent=2)+'\n')
for r in res:
 print(r['gene'],'indels',r['indels'])
 for p in r['base_counts']:
  ref=p['reference'];counts=p['allele_fragment_counts'];alts={a:n for a,n in counts.items() if a!=ref and n>=3}
  if alts:print(p['position_1based'],ref,counts,p['allele_read_strands'])
 print('refdepth range',min(sum(x['allele_fragment_counts'].values()) for x in r['base_counts']),max(sum(x['allele_fragment_counts'].values()) for x in r['base_counts']))
