"""Bounded LATS1/LATS2 candidate discovery, not a clinical somatic caller."""
import pathlib,json,urllib.request,collections,csv,time,hashlib
import pysam
ROOT=pathlib.Path(__file__).resolve().parents[3];OUT=pathlib.Path(__file__).resolve().parent
SRC=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853');OLD=ROOT/'work/oct1-analysis'
ANN=json.loads((OUT/'target-annotation.json').read_text())['genes']
def merge(a):
 out=[]
 for s,e in sorted(a):
  if out and s<=out[-1][1]:out[-1][1]=max(e,out[-1][1])
  else:out.append([s,e])
 return out
refs={}
for gene in ['LATS1','LATS2']:
 g=ANN[gene];f=OUT/f'{gene}.hg38-reference.json';url=f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={g['chrom']};start={g['gene_start0']};end={g['gene_end0']}"
 if not f.exists():
  with urllib.request.urlopen(url,timeout=60) as r:raw=json.load(r)
  assert len(raw['dna'])==g['gene_end0']-g['gene_start0']
  raw['source_url']=url;f.write_text(json.dumps(raw)+'\n')
 refs[gene]=json.loads(f.read_text())
candidates=[];stats=[];start=time.time()
with pysam.AlignmentFile(str(SRC/'DNA_TN26-279853.bam'),'rb',index_filename=str(OLD/'DNA_TN26-279853.bam.bai')) as bam:
 for gene in ['LATS1','LATS2']:
  g=ANN[gene];regions=merge([(s-4,e+4) for s,e in g['selected_CDS']]);ref=refs[gene]['dna'].upper();rs=g['gene_start0'];covered={};indel_stats=[]
  for s,e in regions:
   for col in bam.pileup(g['chrom'],s,e,truncate=True,stepper='nofilter',max_depth=100000,min_base_quality=0,ignore_overlaps=False,compute_baq=False):
    p=col.reference_pos;rb=ref[p-rs];frag=collections.defaultdict(set);fstrands=collections.defaultdict(set);clean=collections.defaultdict(set);inserts=collections.defaultdict(set);indstrands=collections.defaultdict(set)
    for pr in col.pileups:
     r=pr.alignment
     if r.flag & (4|256|512|1024|2048) or r.mapping_quality<20 or pr.is_del or pr.is_refskip:continue
     q=pr.query_position
     if q is None or r.query_qualities is None or r.query_qualities[q]<20:continue
     key=(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name);base=r.query_sequence[q];frag[key].add(base);fstrands[base].add('R' if r.is_reverse else 'F')
     if not any(op==4 for op,n in r.cigartuples) and min(q-r.query_alignment_start,r.query_alignment_end-1-q)>=5:clean[base].add(key)
     if pr.indel:
      d=pr.indel
      if d>0:
       if q+d>=len(r.query_sequence) or min(r.query_qualities[q:q+d+1])<20:continue
       a=(rb,rb+r.query_sequence[q+1:q+d+1])
      else:a=(ref[p-rs:p-rs+1-d],rb)
      inserts[a].add(key);indstrands[a].add('R' if r.is_reverse else 'F')
    fc=collections.Counter(next(iter(a)) for a in frag.values() if len(a)==1);total=sum(fc.values());covered[p]=total
    for b,n in fc.items():
     if b!=rb and b in 'ACGT' and n>=10 and total>=20 and n/total>=.05:
      candidates.append({'gene':gene,'chrom':g['chrom'],'pos1':p+1,'type':'SNV','ref':rb,'alt':b,'alt_queryname_fragments':n,'assessable_anchor_fragments':total,'screen_fraction':n/total,'support_on_both_read_orientations':len(fstrands[b])==2,'clean_alt_fragments_no_softclip_or_end':len(clean[b]),'screen_only':True})
    for (r,a),keys in inserts.items():
     if len(keys)>=10 and total>=20 and len(keys)/total>=.05:
      candidates.append({'gene':gene,'chrom':g['chrom'],'pos1':p+1,'type':'CIGAR_indel_not_normalized','ref':r,'alt':a,'alt_queryname_fragments':len(keys),'assessable_anchor_fragments':total,'screen_fraction':len(keys)/total,'support_on_both_read_orientations':len(indstrands[(r,a)])==2,'clean_alt_fragments_no_softclip_or_end':None,'screen_only':True})
  totalbases=sum(e-s for s,e in regions);ds=[covered.get(p,0) for s,e in regions for p in range(s,e)]
  stats.append({'gene':gene,'selected_transcript':g['selected_transcript'],'CDS_plus4bp_union_bases':totalbases,'minimum_strict_queryname_fragment_depth':min(ds),'bases_ge20_fragment_depth':sum(d>=20 for d in ds),'bases_ge50_fragment_depth':sum(d>=50 for d in ds),'candidates_over_screen_threshold':sum(r['gene']==gene for r in candidates),'regions0':regions})
  print(gene,stats[-1],flush=True)
(OUT/'lats-exploratory-screen.json').write_text(json.dumps({'method':'DNA, MAPQ/BQ>=20, exclude duplicate/QCfail/secondary/supplementary/unmapped; mate collapse by RG+queryname. Examine all bases in MANE CDS +/-4bp. Candidate threshold >=10 alt fragments, >=5% of assessable anchor fragments, depth>=20. Not validated somatic calling; germline, pseudogene, complex indel and structural/splice changes not comprehensively assessed. Indels are raw CIGAR representations and do not use deletion-length-specific reference denominators; candidate screening only.','reference':'Public hg38 UCSC genomic sequence; exact URLs stored in gene reference JSONs.','gene_stats':stats,'candidates':candidates,'elapsed_seconds':time.time()-start},indent=2)+'\n')
if candidates:
 with (OUT/'lats-exploratory-candidates.tsv').open('w') as f:
  w=csv.DictWriter(f,list(candidates[0]),delimiter='\t');w.writeheader();w.writerows(candidates)
