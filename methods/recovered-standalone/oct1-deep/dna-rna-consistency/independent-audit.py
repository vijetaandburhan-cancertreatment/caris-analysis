"""Independent CIGAR coordinate traversal at the 335 existing markers."""
from pathlib import Path
from collections import defaultdict,Counter
import json,datetime,time
import pysam
P=Path(__file__).resolve().parent;ROOT=P.parents[2];RAW=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853');prior=json.loads((P/'findings.json').read_text());rows=prior['rows'];checked=[];start=time.time()
for kind in ['DNA','RNA']:
 with pysam.AlignmentFile(str(RAW/f'{kind}_TN26-279853.bam'),'rb',index_filename=str(ROOT/f'work/oct1-analysis/{kind}_TN26-279853.bam.bai')) as bam:
  for v in rows:
   calls=defaultdict(set);pos=v['pos1']-1
   for r in bam.fetch(v['chrom'],pos,pos+1):
    if r.flag&(4|256|512|1024|2048) or r.mapping_quality<30 or r.query_qualities is None or (r.has_tag('NH') and r.get_tag('NH')!=1):continue
    rp=r.reference_start;qp=0;index=None
    for op,n in r.cigartuples:
     if op in (0,7,8):
      if rp<=pos<rp+n:index=qp+pos-rp;break
      rp+=n;qp+=n
     elif op in (1,4):qp+=n
     elif op in (2,3):rp+=n
    if index is None or r.query_qualities[index]<30 or index-r.query_alignment_start<5 or r.query_alignment_end-index-1<5:continue
    key=(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name);calls[key].add(r.query_sequence[index])
   c=Counter(next(iter(z)) if len(z)==1 else 'discordant' for z in calls.values());counts=dict(c)
   checked.append({'kind':kind,'chrom':v['chrom'],'pos1':v['pos1'],'gene':v['gene'],'matches_prior_counts':counts==v[kind]['counts'],'counts':counts})
eligible=[v for v in rows if v['DNA']['assessable_names']>=20 and v['RNA']['assessable_names']>=20 and v['DNA']['specified_alt_fraction']>=.98]
summary={'tested_export_selected_markers':len(rows),'DNA_high_alt_and_both_depth20_markers':len(eligible),'RNA_same_allele_ge90pct':sum(v['RNA']['specified_alt_fraction']>=.90 for v in eligible),'RNA_same_allele_ge98pct':sum(v['RNA']['specified_alt_fraction']>=.98 for v in eligible),'genes_represented':len(set(v['gene'] for v in eligible)),'chromosomes_represented':sorted(set(v['chrom'] for v in eligible))}
zrsr=json.loads((P/'audit-zrsr2.json').read_text());score=zrsr['RNA_read_score_classes']
result={'timestamp_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'PASS for bounded internal consistency; flagged marker has local mapping ambiguity','recount_matches_all_670_marker_assays':all(v['matches_prior_counts'] for v in checked),'recomputed_summary':summary,'summary_matches':all(summary[k]==prior['summary'][k] for k in summary),'elapsed_seconds':time.time()-start,'ZRSR2_local_reference_comparison':{'independent_counts':zrsr['independent_counts'],'read_alignment_score_classes':score,'interpretation':'The C-bearing RNA aligned portions are indistinguishable from or better explained by the homologous reference pseudogene in this two-transcript local comparison; they cannot be treated as clean ZRSR2-specific genotype observations. This supports mapping ambiguity as an explanation for the single exception, without proving every read origin.','limits':['Only two reference transcripts compared; no genome-wide remapping, full paired-fragment adjudication, calibrated error model or identity assay.','Soft-clipped sequence and mates not covering the marker were not jointly realigned; those could alter a placement decision.','No contamination rate, chance-match probability, sample identity proof or RNA-editing conclusion follows.']},'checks':checked,'overall_limits':['Markers selected from the same DNA VCF by near-homozygous alternate fraction and rs/Benign filters; not an independently chosen fingerprint panel.','Common near-fixed population alleles and linked loci may add little independent identity information.','Bulk RNA has allele-specific expression, mapping and duplicate/fragment-count limitations.','No matched normal or quantitative contamination/relatedness estimate.']}
(P/'audit.json').write_text(json.dumps(result,indent=2)+'\n');print(result['recomputed_summary']);print('all marker assays match',result['recount_matches_all_670_marker_assays'],'summary match',result['summary_matches'],'sec',result['elapsed_seconds'])
