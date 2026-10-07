"""Summarize bounded genomic results with explicit evidence limitations."""
import pathlib,json,csv,collections
B=pathlib.Path(__file__).resolve().parent;ROOT=B.parents[2];OLD=ROOT/'work/oct1-analysis'
x=json.loads((OLD/'variants.xlsx-records.json').read_text())
coverage=list(csv.DictReader((B/'coding-coverage-summary.tsv').open(),delimiter='\t'))
expr={r['Gene']:float(r['TPM']) for r in csv.DictReader(open('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853/RNA_TN26-279853.geneTPM_nodup.csv'))}
genes=['BAP1','MTAP','CDKN2A','CDKN2B','NF2','TP53','SETD2','TERT','LATS1','LATS2','PBRM1','RB1','STAG2','SMARCA2','B2M','JAK1','JAK2','FGFR3','KDM6A','ARID1A','VHL','PAX8']
rows=[]
for gene in genes:
 tests=[r for r in x if r['Biomarker']==gene];cs=[r for r in coverage if r['gene']==gene]
 d={'gene':gene,'Caris_small_variant_result':'; '.join(sorted(set(r['Test Result'] for r in tests if r['Technology']=='Hybrid Exome'))) or 'not explicitly reported','Caris_deletion_result':'; '.join(sorted(set(r['Test Result'] for r in tests if r['Technology']=='Hybrid Exome CND'))) or 'not explicitly reported','Caris_gene_TPM':expr.get(gene),'clinical_negative_warning':'A report negative is assay-limited; depth is not validated sensitivity and does not exclude structural, epigenetic, subclonal or uncaptured changes.'}
 for r in cs:
  label='MANE_or_selected' if r['scope']=='selected_transcript_CDS' else 'all_CDS_union'
  for key in ['bases','median_read_depth','min_read_depth','bases_ge20']:
   source='bases_ge_20' if key=='bases_ge20' else key;d[label+'_'+key]=r[source]
  d[label+'_fraction_ge20']=int(r['bases_ge_20'])/int(r['bases'])
 rows.append(d)
with (B/'gene-evidence-and-coverage.tsv').open('w') as f:
 w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
ann=json.loads((B/'lats-conditional-annotation.json').read_text());read=json.loads((B/'lats-independent-haplotype-recount.json').read_text())
new=[]
for v in ann['candidates']:
 if v['conditional_protein']=='synonymous':continue
 d={k:v[k] for k in ['gene','chrom','pos1','ref','alt','transcript','conditional_coding','conditional_protein']}
 for kind in ['DNA','RNA']:
  r=next(r for r in read['results'] if r['kind']==kind and r['gene']==v['gene'] and r['pos1']==v['pos1'])
  for k in ['alt','ref','other','discordant']:d[kind+'_'+k+'_fragments']=r['fragments'].get(k,0)
  for k in ['alt','ref']:d[kind+'_'+k+'_reads']=r['reads'].get(k,0)
  d[kind+'_alt_clean_supporting_fragments']=r['fragments_with_clean_no_softclip_and_no_end_support'].get('alt',0)
 if v['pos1']==20988693:d.update(public_id='rs2770928',public_gnomADe_observed_allele_frequency=.8806,research_priority='deprioritize as tumor-specific candidate; common population allele')
 elif v['pos1']==20988809:d.update(public_id='rs558614',public_gnomADe_observed_allele_frequency=.7906,research_priority='deprioritize as tumor-specific candidate; common population allele')
 else:d.update(public_id='no exact match returned by limited Ensembl overlap lookup',public_gnomADe_observed_allele_frequency='not established',research_priority='confirm exact variant and normal status; clinical functional/actionability interpretation required')
 d.update(Caris_VCF_or_workbook_explicit_variant='absent',somatic_status='unknown; no matched normal',clinical_classification='not established by this exploratory analysis',protein_status='conditional reference-transcript prediction only',RNA_limit='query-name fragments, not UMI independent molecules; no protein/presentation validation')
 d['separate_implementation_audit']='passed exact DNA/RNA counts and reference-transcript protein mapping; same original data/reference, not clinical validation' if v['pos1'] in [149695191,20975259] else 'not separately audited'
 new.append(d)
with (B/'new-LATS-candidate-evidence.tsv').open('w') as f:
 w=csv.DictWriter(f,list(new[0]),delimiter='\t');w.writeheader();w.writerows(new)
sources=[
 {'purpose':'Reference GTF','url':'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz'},
 {'purpose':'Allele-specific CN, purity/ploidy and copy-neutral LOH methods','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5027494/'},
 {'purpose':'Primary MTAP protein/copy concordance study, heterozygous vs homozygous distinctions','url':'https://pubmed.ncbi.nlm.nih.gov/33387364/'},
 {'purpose':'Independent orthogonal MTAP methods study','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11823465/'},
 {'purpose':'Primary paper documenting exact hg38 TERT hotspot coordinates','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC11914184/'},
 {'purpose':'Primary LATS2 rs558614 germline polymorphism study','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC4792794/'},
 {'purpose':'Official public population annotations rs2770928','url':'https://rest.ensembl.org/variation/human/rs2770928?pops=1;content-type=application/json'},
 {'purpose':'Official public population annotations rs558614','url':'https://rest.ensembl.org/variation/human/rs558614?pops=1;content-type=application/json'},
]
(B/'sources.json').write_text(json.dumps(sources,indent=2)+'\n')
summary={
 'scope':'Bounded deep pass over resident Caris files; research candidates and negative-evidence assessment, not a clinically validated reanalysis or new diagnosis.',
 'MTAP':{'Caris_CND':'MTAP deleted; export CNA value0.80. CDKN2A1.30/CDKN2B0.90 also deleted; values not treated as calibrated integer tumor copy counts.','Tata_Sept18':'Visually verified original final histopathology page1: IHC performed on both brain and ?lymph-node blocks shows p16 and MTAP protein loss in tumor cells. Independent orthogonal support for functional MTAP loss.','gene_coordinate_cluster':'Three genes lie within hg38 chr9:21802636–22009305 (~207kb gene-span envelope), not a measured continuous deletion interval or breakpoint.','depth':'Canonical MTAP coding depth57–305x; median221x. Residual reads/TPM4.877542 do not refute tumor-cell protein loss; normal cells/heterogeneity/capture affect bulk measurements. No new homozygous-loss call.'},
 'BAP1':{'second_hit':'Only reported small variant is R57fs; Caris deletion not detected. Sparse chr3p marker allele imbalance is strong and reproducible, supports reviewing allele-specific copy number/LOH, not a confirmed BAP1 second hit.','allelic_pattern':'Six selected dbSNP SNPs on chr3:0–80Mb have median minor allele fraction0.2185 (range0.189–0.233), vs21 markers100Mb–end median0.385. Markers sparse/ascertainment-biased, no matched normal or physical phasing to BAP1.','why_not_copy_neutral_call':'Both hemizygous loss and copy-neutral LOH with different purity can generate similar allelic fractions. Allelic imbalance alone cannot distinguish copy state, somatic/germline status or clonality.'},
 'new_candidates':new,
 'new_candidate_independent_audit':'Root separate CIGAR-walk/haplotype implementation exactly reproduces both research-priority loci DNA/RNA counts, and separate exon-edit/strand-reversed translation agrees with GENCODE v37 protein FASTA. See lats-root-independent-audit.json/log/py. Same original data/reference; not independent biological or clinical confirmation.',
 'negative_evidence':{'TERT_hotspots':'G>A absent from1295113 (1710G reads,1391G fragments) and1295135 (1647G reads,1295G fragments), strict quality/flag filters. Better-supported narrow negative, not a diagnostic exclusion.','NF2_TP53_SETD2':'MANE CDS fully covered at high strict read depth, compatible with reported small-variant negatives. Alternate-transcript/structural/epigenetic and assay-sensitivity limits remain.','LATS1_LATS2':'No explicit Caris workbook/VCF gene result. Raw-read candidate discovery found two research-priority variants plus common SNPs, demonstrating that the export is not exhaustive. Never translate absent VCF gene into wild type.','RB1':'31/2784 MANE CDS bases <20x (minimum7x). Negative needs this local coverage caveat.','indeterminate':'24 workbook small-variant rows are indeterminate. Public-CDS20x coverage does not override the laboratory designation or recover its validated sensitivity.'},
 'recommended_validation_priorities':['Independent implementation audit + clinical laboratory confirmation of LATS1 c.378dup and LATS2 c.2878G>A, then matched-normal exclusion before tumor-specific/neoantigen claims.','Clarify Caris full-exome output scope and request variant evidence/annotation for these loci. Report omission may be scope/reportability, not an error.','Use existing Tata MTAP IHC plus Caris deletion for trial prescreening; confirm exact protocol-required MTAP biomarker definition with study team.','Obtain calibrated allele-specific segmented CN/purity/ploidy and normal data to evaluate BAP1 second hit; do not infer homozygous deletion from raw exon depth.'],
 'sources':sources,
}
if (B/'reported-missense-public-annotation.tsv').exists() and (B/'reported-missense-clinvar.tsv').exists():
 population=list(csv.DictReader((B/'reported-missense-public-annotation.tsv').open(),delimiter='\t'))
 clinical=list(csv.DictReader((B/'reported-missense-clinvar.tsv').open(),delimiter='\t'))
 combined=[]
 for r in population:
  c=next(c for c in clinical if (c['gene'],c['protein'])==(r['gene'],r['protein']))
  combined.append({**r,**{k:v for k,v in c.items() if k not in r}})
 with (B/'reported-missense-population-clinvar-combined.tsv').open('w') as f:
  w=csv.DictWriter(f,list(combined[0]),delimiter='\t');w.writeheader();w.writerows(combined)
 summary['reported_missense_catalog_review']={
  'count':len(population),
  'exact_alleles_with_gnomAD_population_observations':sum(bool(r['gnomADe:ALL_observed_alt_frequency'] or r['gnomADg:ALL_observed_alt_frequency']) for r in population),
  'no_exact_Ensembl_match':['FANCM p.Y537F','MALT1 p.Q299E','WAS p.P370T'],
  'exact_ClinVar_VCV_records':sum(bool(r['exact_ClinVar_accessions']) for r in clinical),
  'direct_ClinVar_summary':'ALK/AXIN2/PRKAR1A: conflicting germline classifications; CDH23/FANCD2/SDHA/two ZFHX3: uncertain significance; PRDM1: not provided. None of the nine exact matching VCV records has an aggregate somatic oncogenicity or somatic clinical-impact classification.',
  'neoantigen_implication':'Absence of candidate peptides from a reference protein FASTA is not tumor specificity: many rare germline proteins are absent from that reference. Fourteen exact alleles have population observations; all17 require patient-matched normal assessment before calling tumor-specific. The three absent catalog hits are not automatically somatic. Germline benign/pathogenic significance and cancer-driver status are different from peptide immunogenicity.',
  'outputs':['reported-missense-population-clinvar-combined.tsv','reported-missense-public-annotation.json','reported-missense-clinvar.json'],
 }
(B/'findings.json').write_text(json.dumps(summary,indent=2)+'\n')
print('Wrote gene-evidence-and-coverage.tsv, new-LATS-candidate-evidence.tsv, findings.json, sources.json')
