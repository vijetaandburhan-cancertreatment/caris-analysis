"""Independent local-haplotype recount of new LATS screening candidates."""
import json,pathlib,collections,statistics
import pysam
ROOT=pathlib.Path(__file__).resolve().parents[3];OUT=pathlib.Path(__file__).resolve().parent;OLD=ROOT/'work/oct1-analysis'
SRC=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
x=json.loads((OUT/'lats-conditional-annotation.json').read_text())
results=[]
for kind in ['DNA','RNA']:
 with pysam.AlignmentFile(str(SRC/f'{kind}_TN26-279853.bam'),'rb',index_filename=str(OLD/f'{kind}_TN26-279853.bam.bai')) as bam:
  for v in x['candidates']:
   if v['conditional_protein']=='synonymous':continue
   ref=json.loads((OUT/f'{v["gene"]}.hg38-reference.json').read_text());rs=ref['start'];dna=ref['dna'].upper();p=v['pos1']-1
   st=p-8;en=p+len(v['ref'])+8;rh=dna[st-rs:en-rs];ah=rh[:8]+v['alt']+rh[8+len(v['ref']):]
   counts=collections.Counter();strands=collections.Counter();f=collections.defaultdict(set);coords=collections.defaultdict(set);clean=collections.defaultdict(set);mapqs=collections.defaultdict(list)
   for r in bam.fetch(v['chrom'],st,en):
    if r.flag & (4|256|512|1024|2048) or r.mapping_quality<20 or r.query_qualities is None:continue
    aps={rp:q for q,rp in r.get_aligned_pairs(matches_only=True) if rp is not None}
    if st not in aps or en-1 not in aps:continue
    qs=aps[st];qe=aps[en-1]+1
    if min(r.query_qualities[qs:qe])<20:continue
    observed=r.query_sequence[qs:qe];call='ref' if observed==rh else 'alt' if observed==ah else 'other'
    key=(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name);counts[call]+=1;f[key].add(call);strands[call+('_reverse' if r.is_reverse else '_forward')]+=1
    coords[call].add((r.reference_start,r.reference_end,r.cigarstring,r.is_reverse,r.is_read1));mapqs[call].append(r.mapping_quality)
    if not any(op==4 for op,n in r.cigartuples) and min(qs-r.query_alignment_start,r.query_alignment_end-qe)>=5:clean[call].add(key)
   fc=collections.Counter(next(iter(a)) if len(a)==1 else 'discordant' for a in f.values())
   result={**v,'kind':kind,'haplotype_start0':st,'haplotype_end0':en,'reads':dict(counts),'fragments':dict(fc),'strand_reads':dict(strands),'distinct_alignment_coordinate_patterns':{a:len(b) for a,b in coords.items()},'fragments_with_clean_no_softclip_and_no_end_support':{a:len(b) for a,b in clean.items()},'median_MAPQ':{a:statistics.median(b) for a,b in mapqs.items()}}
   results.append(result);print(kind,v['gene'],v['pos1'],counts,fc,flush=True)
(OUT/'lats-independent-haplotype-recount.json').write_text(json.dumps({'method':'Compare continuous observed read sequence between exact public reference endpoints spanning candidate +/-8bp against reference vs alternate haplotypes. Captures equivalent local insertion shifts without relying on CIGAR anchor. All bases BQ>=20, MAPQ>=20; exclude unmapped/secondary/supplementary/QCfail/duplicate. Fragment means RG+queryname, not UMI molecule. Other mismatching haplotypes tracked, discordant mates excluded. Fractions differ from initial pileup due stringent full-haplotype spanning/quality requirements.','results':results},indent=2)+'\n')
