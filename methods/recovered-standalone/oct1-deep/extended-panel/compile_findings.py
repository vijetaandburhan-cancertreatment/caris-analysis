import pathlib,json,csv,collections
B=pathlib.Path(__file__).resolve().parent;ROOT=B.parents[2]
main=json.loads((B/'screen-results.json').read_text());nf1=json.loads((B/'NF1-screen-results.json').read_text());x=json.loads((ROOT/'work/oct1-analysis/variants.xlsx-records.json').read_text())
cov=list(csv.DictReader((B/'coding-coverage-summary.tsv').open(),delimiter='\t'))+list(csv.DictReader((B/'NF1-coding-coverage-summary.tsv').open(),delimiter='\t'))
rows=[]
for r in cov:
 if r['scope']!='selected_transcript_CDS':continue
 gene=r['gene'];tests=[t for t in x if t['Biomarker']==gene]
 row={**r,'Caris_small_variant_result':'; '.join(sorted(set(t['Test Result'] for t in tests if t['Technology']=='Hybrid Exome'))),'Caris_deletion_result':'; '.join(sorted(set(t['Test Result'] for t in tests if t['Technology']=='Hybrid Exome CND'))),'research_screen_protein_altering_candidates':sum(v['gene']==gene and v['consequence'] not in ['synonymous','CDS_flank_or_exon_boundary'] for v in main['candidates']+nf1['candidates'])}
 rows.append(row)
with (B/'all14-gene-coverage-and-lab-results.tsv').open('w') as f:
 w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
summary={
 'scope':'13-gene bounded screen plus separate NF1 addendum. GENCODE v37 MANE CDS +/-4 genomic bp, same resident DNA/RNA Caris BAMs. No new clinical diagnosis or clinically validated variant calls.',
 'genes':[r['gene'] for r in rows],
 'principal_result':'No additional omitted protein-altering variant found in this 14-gene extension at stated screening thresholds. Every threshold-level candidate exactly matches an existing supplied Caris VCF record. This does not change the separate research-priority LATS1/LATS2 findings.',
 'screen_thresholds':'SNP initial threshold >=10 alternate reads, >=5% of quality-filtered anchor read depth, total depth>=20. Indels >=10 alternate query-name fragments and>=5% of anchor read depth; preliminary CIGAR representations, not complete local reassembly. BQ/MAPQ>=20; flagged duplicate/secondary/supplementary/QCfail/unmapped reads removed. No formal limit of detection.',
 'candidate_categories':dict(collections.Counter(v['consequence'] for v in main['candidates']+nf1['candidates'])),
 'common_protein_variants':[
  {'gene':'SETD2','protein':'P1962L','public_id':'rs4082155','gnomADe_observed_alt_frequency':.5618,'ClinVar':'VCV000135202.16 Benign'},
  {'gene':'TP53','protein':'P72R','public_id':'rs1042522','gnomADe_observed_alt_frequency':.7163,'ClinVar':'VCV000012351.85 Benign'},
  {'gene':'ERBB2','protein':'P1170A','public_id':'rs1058808','gnomADe_observed_alt_frequency':.6544,'ClinVar':'VCV000134082.12 Benign'}
 ],
 'flagged_candidates_not_promoted':[
  {'gene':'PTEN','protein':'Y336H','Caris_filter':'RV; CI=VUS','DNA_SNP_alt_reads':84,'failure':'Every alternate SNP read has alternate base at the final aligned position (distance0); none spans +/-8bp. Zero exact high-quality alternate local haplotypes in DNA or RNA.','ClinVar':'VCV000836781.12 Uncertain significance; public clinical annotation does not rescue deficient patient read evidence.'},
  {'gene':'ERBB2','protein':'D904V','Caris_filter':'sb.s;Benign','DNA_SNP_alt_reads':134,'failure':'All134 alternate reads reverse-strand;130 fail full-flank BQ20 and4 fail +/-8bp span. Zero exact high-quality alternate local haplotypes in DNA/RNA.','public_population_note':'Rare allele catalog existence is not evidence that this specimen carries a real variant.'},
  {'gene':'PBRM1','variants':'c.996-5-region insertion/deletion','Caris_filter':'R8.1; insertion VUS, deletion Benign','failure':'Existing repeat-region flagged calls outside CDS; not new splice or driver calls.'}
 ],
 'NF1':{'selected_transcript':'ENST00000358273.9','CDS_bases':8517,'minimum_strict_read_depth':239,'median_strict_read_depth':2094,'raw_screen':'Only two existing synonymous SNPs; no protein-altering candidate above thresholds.','Caris':'Wild Type small-variant result row709; Deletion Not Detected row1137.','clinical_research_implication':'Combined RASA1+NF1 inactivation is not established; this pass provides no positive support for it. Do not transfer efficacy of dual-loss models onto this specimen on RASA1 alone.'},
 'coverage_limit':'RB1 has31/2784 selected-CDS bases<20x, minimum7x. Other selected-CDS regions high-depth; alternate-transcript CDS unions include coverage gaps. Public transcript intervals are not the laboratory capture BED.',
 'negative_result_limits':['Screen sensitivity not clinically calibrated; lower-frequency variants may be missed.','No exhaustive SV, long insertion, intronic, promoter (apart from separate TERT-hotspot review), methylation, splice, fusion or protein-loss exclusion.','No normal sample, capture-normalization reference, or purity/ploidy fit.','High coding depth and absence of a candidate do not exclude a cancer type or epigenetic/structural gene inactivation.'],
 'files':['all14-gene-coverage-and-lab-results.tsv','screen-results.json','screen-candidates.tsv','protein-candidate-haplotype-recount.json','flagged-candidate-failure-modes.json','extended-protein-candidate-public-annotation.json','extended-protein-candidate-clinvar.json','NF1-screen-results.json'],
}
(B/'findings.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({'genes':len(rows),'candidate_categories':summary['candidate_categories'],'principal_result':summary['principal_result']},indent=2))
