from pathlib import Path
import json,csv,hashlib,shutil,datetime
R=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct4-hla-allele-support');W=Path('/Users/burhanazeem/Documents/Codex/2026-09-05/finances-plugin-finances-openai-curated-remote-3/work/oct4-followup/hla-allele-support');W.mkdir(parents=True,exist_ok=True)
j=json.load(open(R/'common-patient-counts.json'));a=json.load(open(R/'independent-recount91.json'));mask=json.load(open(R/'nonB-genomic-mask.json'));control=json.load(open(R/'common-public-controls.json'))
assert a['matches_primary_counts'] and a['matches_primary_sequence_families'] and not any(m['k']>=51 for m in mask['new_nonB_matches'])
TS=['B*50:01','B*52:01'];rows=[]
for k in (51,71,91,111,131):
 c=j['totals']['Q30_end5'][str(k)]['matched_positions']
 for t in TS:rows.append({'marker_length_nt':k,'allele_compatible_with_prior_genotype':t,'base_quality_min':30,'read_end_distance_min_nt':5,'matched_CDS_start_positions':control['matched_start_counts'][str(k)],'positive_query_name_pairs':c.get(t,0),'distinct_full_pair_sequences':j['distinct_full_pair_sequences_matched_Q30_end5'][t][str(k)],'cross_allele_conflict_pairs':c.get(';'.join(TS),0),'interpretation':'Conditional bulk RNA support; not unique HLA retyping, tumor retention, independent molecules or calibrated allele abundance'})
with open(R/'matched-window-support.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
co=j['marker_pair_counts']['Q30_end5'];foot=[]
for i,m in enumerate(j['markers']):
 if m['k']==91 and m['matched_position']:
  foot.append({'allele':m['allele'],'CDS_start_1based':m['positions'][0][1]+1,'CDS_end_1based':m['positions'][0][1]+91,'marker_length':91,'pair_count_Q30_end5':co.get(str(i),0),'number_compatible_other_HLA_B_groups':len(m['compatible_two_field'])-1})
with open(R/'marker-footprint91.tsv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(foot[0]),delimiter='\t');w.writeheader();w.writerows(foot)
text='''Conditional HLA-B*50:01 RNA support — 5 October 2026

There is substantial HLA-B*50:01-compatible RNA sequence in the supplied bulk specimen. Given the existing B*50:01/B*52:01 genotype interpretation, this supports continuing the BAP1/HLA-B*50:01 research candidate. It does not prove that the malignant cells retain or express that allele, or that they display the BAP1 peptide.

The primary follow-up used all23,209,264 original RNA read pairs, without genomic alignment or chr6 selection. At344 matched91-nucleotide coding-window starts, requiring every marker base Q>=30 and at least5 nucleotides between the marker and either read end, there were2,132 B*50:01-compatible pairs and1,655 B*52:01-compatible pairs. They contained1,588 and1,299 distinct full read-pair sequences, respectively. Every one of the344 matched windows was observed for each allele; window starts span coding positions1–756 (marker ends through846). These overlapping windows and read-sequence patterns are not independent molecules. A separate four-line FASTQ parser and string-search implementation exactly reproduced the primary91-mer counts and sequence-family counts.

Sensitivity checks retain both patterns: with the same stringent Q30/end5 filters, matched71-mers give2,564/1,944 pairs,111-mers1,875/1,450, and131-mers1,291/1,005, in B50/B52 order. Q25/end5 matched91-mers give2,546/1,933. Counts are conditional sequence-support measurements, not a calibrated expression ratio. In particular, naive short31-mer markers gave misleadingly unequal detection opportunities (only18 gene-specific common-B50 starts versus116 common-B52 starts before matched-window restriction). Matching the same coding starts removes that simple design imbalance, but does not eliminate transcript coverage, quality, amplification, alternate-isoform or genotype uncertainty.

Specificity controls used all32,330 available coding-sequence records in IPD-IMGT/HLA3.46.0, including partial alleles, all234,485 GENCODE37 transcripts, and all23,995 non-HLA-B genomic records from the same HLA catalog (including pseudogenes/records without CDS). No51–131mer main-analysis marker acquired an extra non-B genomic hit; one preliminary31mer did and is excluded from final interpretation. Exact-match, overlapping-match and reverse-complement public controls passed. The two common coding representatives align without gaps over1,089 nucleotides;344 shared91-mer starts have equal detection opportunity in ideal equal-start public fragment simulations across75/100/125/150nt read lengths and200/250/300nt fragment lengths.

Crucial limit: no individual common-reference31–131nt marker, nor the observed within-fragment compatible-group intersections, uniquely excludes every other catalogued HLA-B allele. These sequence blocks are shared with other HLA-B types. Therefore this is conditional support for the previously inferred B50/B52 pair, not independent clinical HLA typing or high-field genotype confirmation. The catalog is the available3.46.0 release, not all later or unknown alleles, and missing sections of partial alleles cannot be excluded by absence from the reference.

Biological limit: bulk tumor RNA includes nonmalignant immune/stromal cells. This analysis cannot assign the B50-compatible transcripts to malignant cells, establish HLA copy retention/LOH, surface expression, antigen processing, peptide presentation, T-cell recognition or vaccine efficacy. The BAP1 mutant transcript and HLA transcript have not been linked to the same individual cell. No matched normal or new tissue is present in this analysis. The strongest justified claim is 'B*50:01-compatible bulk RNA is readily detectable under the existing genotype interpretation.'

Reproducibility: prepare_common_haplotypes.py, count_common_haplotypes.py, mask_nonB_genomic_hla.py, independent_recount91.py; reference and public controls, per-marker counts and local candidate FASTQs retained. Patient raw FASTQs were not modified. Full scan377.66seconds, peak0.86GB RSS; the earlier verbose reference-masking pass peaked1.98GB, and the independent recount peaked1.25GB. This HLA inference/recount ran locally; it made no external patient-data submissions and used no paid compute. Public package pyahocorasick2.2.0 was installed only in this task-local dependency folder. Original sparse31-mer pilot files are historical diagnostics, not the final allele-expression result.

Primary source provenance:
https://github.com/ANHIG/IMGTHLA/tree/3460
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz
https://github.com/RabadanLab/arcasHLA
https://github.com/WojciechMula/pyahocorasick
'''
# Improve visible prose spacing without altering identifiers or exact values.
import re
text=re.sub(r'(?<=[a-zA-Z])(?=\d)', ' ',text)
# Keep formal identifiers and version spellings readable; URLs untouched from original below.
text=text.replace('B*','B*').replace('Q>=','Q>=')
# Restore URLs altered by the generic readability step.
text=text.replace('tree/3460','tree/3460').replace('gencode.v 37','gencode.v37').replace('release_37','release_37')
for old,new in [('BAP 1','BAP1'),('B 50','B50'),('B 52','B52'),('chr 6','chr6'),('Q 30','Q30'),('Q 25','Q25'),('end 5','end5'),('recount 91.py','recount91.py'),('pairs,111','pairs, 111'),('nucleotides;344','nucleotides; 344'),('131mer','131-mer'),('31mer','31-mer')]:text=text.replace(old,new)
text=re.sub(r'(?<=\d)(?=(?:nt|GB|seconds)\b)',' ',text)
(R/'findings.txt').write_text(text)
o={'status':'complete','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'main_result':rows[4:6],'all_sensitivity_rows':rows,'independent_recount_pass':True,'full_original_RNA_pairs':j['counts']['pairs'],'reference_global_specificity':'conditional; other HLA-B allele compatibility unresolved','nonB_genomic_mask':{k:v for k,v in mask.items() if k!='new_nonB_matches'},'tumor_retention_or_LOH_claim':False,'external_patient_submission_by_this_HLA_workstream':False,'paid_compute_cost_usd':0,'main_findings_path':str(R/'findings.txt')};(R/'findings.json').write_text(json.dumps(o,indent=2))
files=['findings.txt','findings.json','matched-window-support.tsv','marker-footprint91.tsv','common-haplotype-summary.json','common-public-controls.json','representative-coordinate-alignment.json','nonB-genomic-mask.json','independent-recount91.json','prepare_common_haplotypes.py','count_common_haplotypes.py','mask_nonB_genomic_hla.py','independent_recount91.py','compile_findings.py','patient-scan.log','common-patient-scan.log','mask-and-control-summary.json','intermediate-compression.json']
manifest=[]
for name in files:
 p=R/name
 if p.exists():
  shutil.copy2(p,W/name);manifest.append({'file':name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(R/'handoff-manifest.json').write_text(json.dumps({'resident_root':str(R),'workspace_handoff_root':str(W),'files':manifest,'large_evidence_retained_resident':['common-patient-counts.json','common-haplotype-markers.json.gz','common-positive.R1.fastq.gz','common-positive.R2.fastq.gz','common-positive-evidence.jsonl.gz','markers-unmasked.json.gz','markers.json','gencode.v37.transcripts.fa.gz']},indent=2));shutil.copy2(R/'handoff-manifest.json',W/'handoff-manifest.json');print(json.dumps(o),flush=True)
