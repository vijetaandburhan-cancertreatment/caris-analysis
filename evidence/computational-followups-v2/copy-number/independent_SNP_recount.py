"""Separate direct CIGAR-position recount, with both-read-end sensitivity."""
import pathlib,csv,json,collections,random,time
import pysam
OUT=pathlib.Path(__file__).resolve().parent;ROOT=OUT.parents[2];SRC=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853');rng=random.Random(20261004)
rows=list(csv.DictReader((OUT/'public-common-SNP-alleles.annotated.tsv').open(),delimiter='\t'));stric=[r for r in rows if r['candidate_heterozygous_stringent']=='1'];selected={}
for arm in sorted(set(r['arm'] for r in stric)):
 a=[r for r in stric if r['arm']==arm]
 for r in rng.sample(a,min(25,len(a))):selected[(r['chrom'],r['pos1'])]=r
for r in stric:
 if (r['chrom']=='chr3' and 51000000<=int(r['pos1'])<=53500000) or (r['chrom']=='chr9' and 19000000<=int(r['pos1'])<=25000000):selected[(r['chrom'],r['pos1'])]=r
bam=pysam.AlignmentFile(str(SRC/'DNA_TN26-279853.bam'),'rb',index_filename=str(ROOT/'work/oct1-analysis/DNA_TN26-279853.bam.bai'));results=[]
for _,r in sorted(selected.items()):
 pos=int(r['pos1'])-1;d=collections.defaultdict(set);clean=collections.defaultdict(set);mq60=collections.defaultdict(set)
 for a in bam.fetch(r['chrom'],pos,pos+1):
  if a.flag&0xF0C or not a.is_proper_pair or a.mapping_quality<30 or a.query_qualities is None:continue
  mp={p:q for q,p in a.get_aligned_pairs(matches_only=True)}
  if pos not in mp:continue
  q=mp[pos]
  if a.query_qualities[q]<25:continue
  b=a.query_sequence[q]
  if b not in 'ACGT':continue
  key=a.query_name;d[key].add(b)
  if a.mapping_quality==60:mq60[key].add(b)
  if not any(op==4 for op,n in a.cigartuples) and min(q-a.query_alignment_start,a.query_alignment_end-1-q)>=5:clean[key].add(b)
 def tally(groups,master=None):
  c=collections.Counter()
  for k,v in groups.items():
   if master is not None and len(master[k])!=1:continue
   if len(v)==1:c[next(iter(v))]+=1
   else:c['discordant']+=1
  return c
 base=tally(d);cl=tally(clean,d);hi=tally(mq60,d)
 res={'chrom':r['chrom'],'pos1':int(r['pos1']),'ref':r['ref'],'alt':r['alt'],'arm':r['arm'],'original_ref':int(r['fragment_ref']),'original_alt':int(r['fragment_alt']),'independent_ref':base[r['ref']],'independent_alt':base[r['alt']],'independent_discordant':base['discordant'],'clean_ref':cl[r['ref']],'clean_alt':cl[r['alt']],'MAPQ60_ref':hi[r['ref']],'MAPQ60_alt':hi[r['alt']]};res['exact_refalt_match']=res['original_ref']==res['independent_ref'] and res['original_alt']==res['independent_alt'];results.append(res)
with (OUT/'independent-SNP-recount.tsv').open('w') as f:w=csv.DictWriter(f,list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
summary={'loci':len(results),'exact_matches':sum(r['exact_refalt_match'] for r in results),'mismatches':[r for r in results if not r['exact_refalt_match']],'selection':'fixed-seed25stringentloci perarm plusallstringentmarkersBAP1regional51–53.5Mb and9p21regional19–25Mb','independence':'separate direct CIGAR position counting implementation, same BAM data/reference; not independentclinical validation','clean_sensitivity':'no softclips, originalbaseBQ>=25,MAPQ>=30,base>=5ntfrombothalignedreadends; discordantmatesexcluded','caveat':'Clean counts are a subset and may have lowerdepth; minorallele ratios rather than totalcount equality shouldbe compared.'};(OUT/'independent-SNP-recount.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
