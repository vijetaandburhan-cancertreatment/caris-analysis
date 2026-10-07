"""Compile evidence without promoting descriptive profiles to CN calls."""
import pathlib,json,csv,statistics,hashlib,datetime,collections
import numpy as np
OUT=pathlib.Path(__file__).resolve().parent;ROOT=OUT.parents[2]
def rd(name):return list(csv.DictReader((OUT/name).open(),delimiter='\t'))
baf=rd('regional-allelic-imbalance-summary.tsv');baf={r['region']:r for r in baf if r['filter']=='candidate_heterozygous_stringent'};depth=rd('chromosome-arm-density-summary.tsv');genes={r['gene']:r for r in rd('selected-CDS-gene-depth.tsv')};snp=rd('public-common-SNP-alleles.annotated.tsv');phase=json.loads((OUT/'BAP1-SNP-frameshift-fragment-phase.json').read_text());audit=json.loads((OUT/'independent-SNP-recount.json').read_text())
counts={'independent_public_SNP_positions_selected':sum(v['selected'] for v in json.loads((OUT/'public-panel-selection.json').read_text())['counts'].values()),'positions_with_pileup_output':len(snp),'base_candidate_heterozygotes':sum(r['candidate_heterozygous_base']=='1' for r in snp),'stringent_candidate_heterozygotes':sum(r['candidate_heterozygous_stringent']=='1' for r in snp),'coding_genes':len(genes),'coding_intervals':len(rd('selected-CDS-interval-depth.tsv')),'direct_CIGAR_audit_loci':audit['loci'],'direct_CIGAR_exact_matches':audit['exact_matches']}
relative=[]
for bins in ['100','250']:
 for pad,mq in [('100','30'),('500','30'),('2000','30'),('500','60')]:
  d={r['arm']:r for r in depth if r['bin_kb']==bins and r['mask_flank']==pad and r['MAPQ_min']==mq}
  relative.append({'bin_kb':int(bins),'mask_flank_bp':int(pad),'MAPQ_min':int(mq),'chr3p_over_chr3q_relative_fragment_density':2**(float(d['chr3p']['adjusted_log2_density_median'])-float(d['chr3q']['adjusted_log2_density_median']))})
direct=rd('independent-SNP-recount.tsv');delta=[]
for r in direct:
 n=int(r['clean_ref'])+int(r['clean_alt']);nr=int(r['original_ref'])+int(r['original_alt'])
 if n>=30:delta.append(abs(min(int(r['clean_ref']),int(r['clean_alt']))/n-min(int(r['original_ref']),int(r['original_alt']))/nr))
source=json.loads((ROOT/'work/oct1-deep/genomics/resident-revalidation.json').read_text())
srcprev=next(r for r in source['files'] if r['name']=='DNA_TN26-279853.bam')
bam=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853/DNA_TN26-279853.bam')
genome=json.loads((OUT/'genome-fragment-profile.method.json').read_text());fragcount=0;records=0
for p in (OUT/'genome-bins').glob('chr*.json'):
 d=json.loads(p.read_text())['counts'];fragcount+=d['accepted_read1_fragments'];records+=d['alignment_records_seen']
findings={
 'created_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'scope':'Expanded research analysis of existing local Caris DNA BAM and analytical export, with public hg38/GENCODE37/1000G reference annotations. No new patient upload or clinical diagnosis. Original sources unchanged.',
 'counts':{**counts,'primary_chromosome_alignment_records_scanned':records,'accepted_read1_fragment_records_for_density':fragcount},
 'BAP1':{
  'finding':'The frameshift is physically linked to the locally overrepresented C allele of a common population SNP. Broad3p and within-gene allelic imbalance are robust. This strengthens a hypothesis of allele-specific alteration involving the mutant haplotype; it is not a confirmed second-hit/LOH/absolute-copy result.',
  'common_SNP':{'hg38':'chr3:52408338 G>C','public_global_AF':.16,'fragment_G':30,'fragment_C':109,'clean_G':26,'clean_C':98,'all_refalt_MAPQ60':True},
  'phase_10bp':next(r for r in phase['results'] if r['flank_bases']==10),
  'regional_allelic_evidence':{k:baf[k] for k in ['chr3p','chr3q','BAP1_plusminus0bp','BAP1_plusminus1000000bp','BAP1_plusminus5000000bp']},
  'relative_density_sensitivity':relative,
  'interpretation_boundaries':['Twenty-one C+deletion pairs do not mean every C haplotype is mutated: seven C+reference and nine G+reference pairs also exist.','All37 linkage observations are paired fragments, not a single read spanning both sites.','Haplotype linkage and tumor-only imbalance cannot distinguish hemizygous loss, copy-neutral LOH, allele-specific gain, mixtures/subclones, or exactsomatic/germline status without more information.','No confirmed biallelic inactivation, protein loss, or new treatment eligibility is claimed.']},
 '9p21':{
  'analytical_export':'Original workbook CND test calls MTAP Deleted/CNA_Value0.80, CDKN2A Deleted/1.30, CDKN2B Deleted/0.90. Those values are retained as vendor outputs, not reinterpreted as integer tumor copies.',
  'local_SNPs':[{'hg38':'chr9:21968160 G>A','reference_fragments':184,'alternate_fragments':198,'clean_ref_alt':[172,186]},{'hg38':'chr9:21993965 T>C','reference_fragments':75,'alternate_fragments':95,'clean_ref_alt':[67,85]}],
  'finding':'The two informative CDKN2A SNPs are approximately balanced despite substantial imbalance elsewhere on9p/9q. Such bulk read evidence can be compatible with normal-cell contribution to a tumor-deleted locus, but it is not specific for homozygous tumor loss.',
  'depth_trough':'An exploratory broad off-exon density trough is visible around22.1–30Mb onchr9, largely downstream of the21.8–22.0Mb MTAP/CDKN2A/B cluster. It must not be presented as the cluster deletion interval or exact breakpoints.',
  'not_established':['Calibrated absolute copy number, homogeneous homozygous loss, purity-adjusted deletion extent and exact breakpoints.','Regional discordant/split alignment inspection did not reveal a supported MTAP/CDKN2A/B-spanning deletion junction; exome capture means this does not exclude one.']},
 'RASA1':{'gene_depth':genes['RASA1'],'finding':'Chromosome5 contains multiple allelic-balance/depth patterns. There are no qualifying common-SNP heterozygotes inside RASA1 or within1Mb under these filters, so local LOH/second-hit status remains unestablished. A distant5Mb-window summary must not substitute for a gene-specific call.','VAF_limit':'The expressed c.747delG truncating allele remains supported by prior read analysis; VAF alone does not establish somatic status, biallelic loss or clonality.'},
 'independent_phase_audit': json.loads((OUT/'BAP1-SNP-phase-independent-audit.json').read_text()) if (OUT/'BAP1-SNP-phase-independent-audit.json').exists() else None,
 'gene_depth_examples':{g:genes[g] for g in ['BAP1','RASA1','MTAP','CDKN2A','CDKN2B','APC']},
 'negative_control_and_sensitivity':{'separate_recount':audit,'absolute_MAF_change_after_clean_both_read_end_filter':{'median':statistics.median(delta),'p95':float(np.quantile(delta,.95)),'max':max(delta),'loci_with_clean_depth_ge30':len(delta)},'synthetic_control':'Passed original-quality/MAPQ/duplicate/supplementary filters, paired-name collapse, conflictingmate exclusion, and zero-depth parsing. Execution validation, not clinical sensitivity/specificity.'},
 'missing_calibration':['Matched normal DNA','Exact private Kapa capture target BED','Vendor reference-normal/bias model','Original CN segments/log2 ratios/BAF/purity/ploidy model','Exact private hg38.pU2AF1_Y_PAR_masked.noHap.fa; primaryhg38 coordinates used here'],
 'next_most_decisive_data':['Caris allele-specific CN/LOH segmentation and purity/ploidy atBAP1/3p, MTAP-CDKN2A/B/9p21 andRASA1/5q.','Caris calibrated MTAP absolute-copy interpretation/thresholds, especially whether homozygousloss is supported and which assayconfirmed it.','Clinician-directed matched-normal testing plus orthogonal BAP1 variant/copy/LOH evaluation when clinically appropriate; this cannot be replaced by the present tumor-only model.'],
 'source_integrity':{'BAM_path':str(bam),'bytes_now':bam.stat().st_size,'previous_verified_SHA256':srcprev['previous_sha256'],'current_size_matches_verified':bam.stat().st_size==srcprev['bytes_now'],'fresh_full_BAM_hash':'not repeated; complete primarychromosome readstream plusindexedrecounts succeeded','public_references':json.loads((OUT/'public/manifest.json').read_text()),'1000G_published_MD5':json.loads((OUT/'public/1000G-download-verification.json').read_text())},
 'primary_method_sources':['https://www.internationalgenome.org/announcements/Variant-calls-from-1000-Genomes-Project-data-on-the-GRCh38-reference-assemlby/','https://www.gencodegenes.org/human/release_37.html','https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/','https://cnvkit.readthedocs.io/en/stable/baf.html','https://cnvkit.readthedocs.io/en/stable/pipeline.html'],
 'resource_use':{'peak_RSS_genome_profile_bytes':max(json.loads(p.read_text())['maxrss_bytes'] for p in (OUT/'genome-bins').glob('chr*.json')),'genome_profile_seconds':genome['elapsed_seconds'],'allele_recount_peakRSS_bytes':json.loads((OUT/'common-SNP-recount.method.json').read_text())['maxrss_bytes'],'public_download_bytes':sum(r['bytes'] for r in json.loads((OUT/'public/manifest.json').read_text()))}
}
# Keep the concise machine summary free of patient query-name hashes/details;
# full phase evidence remains in its dedicated JSON for technical review.
findings['BAP1']['phase_10bp']={k:v for k,v in findings['BAP1']['phase_10bp'].items() if k!='joint_fragment_details'}
if findings['independent_phase_audit']:
 findings['independent_phase_audit']={k:v for k,v in findings['independent_phase_audit'].items() if k!='results'}
(OUT/'findings.json').write_text(json.dumps(findings,indent=2)+'\n')
text=f'''Expanded tumor-only copy-number / allele-balance research review — 4 October 2026

Most useful new result
The BAP1 frameshift is linked to the locally overrepresented SNP haplotype. At hg38 chr3:52408338 G>C, there are109 C versus30 G query-name fragments (clean98 versus26; all originalref/alt observations also passMAPQ60). This SNP lies212bp from the frameshift anchor. With10bp of exact high-quality reference sequence on each side of the deletion,37 informative proper-pair fragments yield:
  C + BAP1deletion: 21 (clean19)
  C + reference at the deletion locus: 7 (clean4)
  G + reference at the deletion locus: 9 (clean8)
  G + BAP1deletion: 0
All37 represent distinct paired-alignment patterns. They are paired-fragment observations, not single reads spanning both loci. The C+reference observations matter: not every C haplotype is mutated. The result supports local cis linkage to the overrepresented allele and a stronger allele-specific alteration hypothesis; it does not establish somatic status, a lost wild-type allele, or clinically confirmed biallelic BAP1inactivation.

The denser regional evidence
The independent public panel contains {counts['independent_public_SNP_positions_selected']:,} predefined common-SNP positions across six chromosomes, with {counts['positions_with_pileup_output']:,} positions producing pileup output and {counts['stringent_candidate_heterozygotes']:,} stringent candidate heterozygotes. There are439 qualifying sites on3p (median minorAF 0.2184), compared with582 on3q (0.3750). BAP1±1Mb has57 markers, median 0.2242. These observations extend the prior six-marker result substantially. No matched normal confirms constitutional heterozygosity.
Relative off-exon DNA fragment density on3p versus3q is0.746–0.764 in100kb bins across exon-exclusion flanks100/500/2000bp. This is a relative within-sample density comparison, not an absolute copy ratio to diploid normal. A region with one retained copy, copy-neutral LOH, allele-specific gain or a mixture can share a similar minorAF at different purity/ploidy. The explicit nonidentifiability table preserves those alternatives instead of selecting a copy state.

MTAP / CDKN2A/B
The original analytical workbook's deletion calls remain the strongest supplied vendor copy-number evidence: MTAP Deleted/CNA_Value0.80; CDKN2A Deleted/1.30; CDKN2B Deleted/0.90. We do not turn these noninteger values into tumor copy counts. Two local CDKN2A SNPs are nearly balanced (184/198 and75/95 reference/alternate fragments), unlike many surrounding chromosome9markers. Normal-cell reads at a tumor-deleted locus are one explanation, but this is not proof of homogeneous homozygous tumor loss. A broader density trough near22.1–30Mb is largely DOWNSTREAM of the21.8–22.0Mb gene cluster and is not an inferred deletion boundary. The bounded discordant/split-read screen did not establish a deletion junction spanning these genes; capture gaps limit this negative observation.

RASA1
The earlier expressed frameshift remains supported. The expanded scan finds no qualifying common-SNP heterozygote within the gene or±1Mb, so local second-hit/LOH status is still unresolved. We measured all3,141 selected coding bases at>=20x, with median282x and mean390.3x. This supports sequence assessability, not a copy-number classification. More distant5q imbalance cannot be assigned to RASA1without additional local evidence.

Why exon depth is not an absolute answer
The expanded coding review covers818genes and11,181selected-CDS intervals. Mean base depth ranges from202.9x atMTAP and390.3x atRASA1 to4447.4x atBAP1; CDKN2A is1845.5x andCDKN2B986.9x. The metadata identify an exome-plus-720-gene workflow, and its private target BED and normal/bias model are absent. Even genes labeled clinical in the export have widely different depth. Within-sample GC/repeat adjustment cannot correct an unknown locus-specific capture baseline. We therefore do not produce integer CN, calibrated homozygous-loss or purity estimates from these depths.

Validation and limits
A separate analyst also reproduced all BAP1 phase counts and fragment identities, verified reciprocal mate coordinates and TLEN, and found no pair-consistency failures. A separate direct CIGAR-position implementation exactly reproduced all388audited SNP reference/alternate fragment counts. Removing soft-clipped reads and bases within5ntof either aligned read end changed minorAF by a median0.0067 across the retained audit sites. Synthetic filter/parser tests passed. Public1000Gdownload matches the publisher's MD5. These are analysis-execution checks, not a clinically validated assay. Fragments are grouped by query name, not UMI-defined independent molecules. No patient data were sent externally, no originals changed, and no clinical source excluded by the user's preference was used.

Most useful next evidence
Request the vendor's allele-specificCN/LOH segments, BAF/log2 ratios, purity/ploidy and target/reference-calibration information for3p/BAP1,9p21/MTAP-CDKN2A/B and5q/RASA1. Ask whether its MTAPresult meets a validated homozygousloss definition, rather than inferring that from0.80. Clinician-directed matched-normal/orthogonal testing is the decisive way to clarify constitutional status and BAP1second-hit interpretation.

Reproducibility
Allscripts, per-position/per-interval tables, count controls, methods and publicreferenceprovenance are in this directory. Main figures: regional-common-SNP-allele-balance.png and regional-off-exon-relative-density.png. Detailed findings.json records exact evidence and missing calibration. BAP1-SNP-frameshift-fragment-phase.json preserves the physical-fragment support; no individual reads are reproduced in this narrative.
'''
# Human-readable spacing, including measurements next to gene names, matters.
replacements={'there are109':'there are 109','versus30':'versus 30','clean98':'clean 98','versus26':'versus 26','passMAPQ60':'pass MAPQ60','originalref/alt':'original ref/alt','lies212bp':'lies 212 bp','With10bp':'With 10 bp','deletion,37':'deletion, 37','BAP1deletion':'BAP1 deletion','clean19':'clean 19','clean4':'clean 4','clean8':'clean 8','All37':'All 37','BAP1inactivation':'BAP1 inactivation','There are439':'There are 439','on3p':'on 3p','with582':'with 582','on3q':'on 3q','minorAF':'minor allele fraction','BAP1±1Mb':'BAP1 ±1 Mb','has57':'has 57','six-marker':'six-marker','3p versus3q':'3p versus 3q','is0.746':'is 0.746','in100kb':'in 100 kb','flanks100/500/2000bp':'flanks 100/500/2000 bp','CNA_Value0.80':'CNA_Value 0.80','markers':'markers','(184/198 and75/95':'(184/198 and 75/95','chromosome9markers':'chromosome 9 markers','near22.1–30Mb':'near 22.1–30 Mb','of21.8–22.0Mb':'of 21.8–22.0 Mb','the21.8–22.0Mb':'the 21.8–22.0 Mb','or±1Mb':'or ±1 Mb','all3,141':'all 3,141','at>=20x':'at ≥20×','median282x':'median 282×','mean390.3x':'mean 390.3×','distant5q':'distant 5q','RASA1without':'RASA1 without','covers818genes and11,181selected-CDS':'covers 818 genes and 11,181 selected-CDS','from202.9x atMTAP and390.3x atRASA1 to4447.4x atBAP1':'from 202.9× at MTAP and 390.3× at RASA1 to 4447.4× at BAP1','is1845.5x andCDKN2B986.9x':'is 1845.5× and CDKN2B 986.9×','all388audited':'all 388 audited','within5ntof':'within 5 nt of','median0.0067':'median 0.0067','Public1000Gdownload':'Public 1000G download','allele-specificCN/LOH':'allele-specific CN/LOH','for3p/BAP1,9p21/MTAP-CDKN2A/B and5q/RASA1':'for 3p/BAP1, 9p21/MTAP-CDKN2A/B and 5q/RASA1','MTAPresult':'MTAP result','homozygousloss':'homozygous loss','from0.80':'from 0.80','BAP1second-hit':'BAP1 second-hit','Allscripts':'All scripts','publicreferenceprovenance':'public reference provenance'}
for a,b in replacements.items():text=text.replace(a,b)
(OUT/'findings.txt').write_text(text)
print(json.dumps(counts,indent=2))
