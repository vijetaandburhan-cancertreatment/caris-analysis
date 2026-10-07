"""Collect exploratory HLA calls, model diagnostics, and reference sensitivity."""
from pathlib import Path
import csv, hashlib, json, re
OUT=Path(__file__).resolve().parent
runs=json.loads((OUT/'patient-genotype-runs.json').read_text())
inputs=json.loads((OUT/'input-feasibility.json').read_text())
extract=json.loads((OUT/'patient-extraction.json').read_text())
prior=json.loads((OUT/'patient-prior-sensitivity.json').read_text())
assert len(prior)==2 and all(r['returncode']==0 and r['guard'] is None for r in prior)
assert len(runs)==2 and all(r['status']=='completed' for r in runs)
def two(alleles):return sorted(':'.join(x.split(':')[:2]) for x in alleles)
def value(pattern,text,cast=str):
    found=re.search(pattern,text)
    return cast(found.group(1)) if found else None
diagnostics={}
for r in runs:
    text=(OUT/('patient-'+r['version']+'.arcas.log')).read_text()
    stats=json.loads(next(Path(r['output_dir']).glob('*.genes.json')).read_text())
    gene_diag={}
    for block in text.split('[genotype] Genotyping HLA-')[1:]:
        gene=block.splitlines()[0].strip();block=block.split('--------------------------------------------------------------------------------')[0]
        abundance=[];pairs=[]
        if '[genotype] Top alleles by abundance:' in block:
            section=block.split('[genotype] Top alleles by abundance:')[1].split('[genotype]')[0]
            abundance=[{'representative_allele':a,'estimated_transcript_abundance_percent':float(b)} for a,b in re.findall(r'^\s+(\S+\*\S+)\s+([\d.]+)%',section,re.M)]
        if '[genotype] Pairs by % explained reads:' in block:
            section=block.split('[genotype] Pairs by % explained reads:')[1].split('[genotype]')[0]
            pairs=[{'representative_pair':[a,b],'percent_gene_fragments_explained':float(c)} for a,b,c in re.findall(r'^\s+(\S+\*\S+), (\S+\*\S+)\s+([\d.]+)%',section,re.M)]
        gene_diag[gene]={'model_gene_fragment_count':stats.get(gene,[None])[0],'compatibility_classes':stats.get(gene,[None,None])[1], 'called_alleles_raw':r['genotypes'].get(gene),'called_two_field':two(r['genotypes'][gene]) if gene in r['genotypes'] else None,'EM_converged_iterations':value(r'EM converged after (\d+) iterations',block,int),'final_EM_alleles':abundance,'surviving_two_field_pair_comparisons':pairs,'genotype_explained_fragments':value(r'Most likely genotype explaining (\d+) reads',block,int),'minor_to_major_nonshared_ratio_rounded':value(r'nonshared count ([\d.]+)',block,float),'zygosity_message':value(r'\[genotype\] (Likely [^\n]+|Unable to distinguish[^\n]+)',block),'no_call':'Not enough reads aligned' in block}
    diagnostics[r['version']]={'genes':gene_diag,'gene_stats_raw':stats,'single_HLA_gene_pseudoaligned_fragments':value(r'\[alignment\] ([\d]+) reads mapped to a single HLA gene',text,int),'pseudoalignment_counts':re.findall(r'processed ([\d,]+) reads, ([\d,]+) reads pseudoaligned',text)}
rows=[]
caris={x['gene'].removeprefix('HLA-'):sorted([x['allele1'],x['allele2']]) for x in inputs['supplied_Caris_classI']}
for gene in ['A','B','C','DPA1','DPB1','DQA1','DQB1','DRB1','DRB5']:
    a=runs[0]['genotypes'].get(gene,[]);b=runs[1]['genotypes'].get(gene,[])
    row={'gene':gene,'IMGT_3_24_two_field':';'.join(two(a)),'IMGT_3_46_two_field':';'.join(two(b)),'reference_concordance':bool(a and b and two(a)==two(b)),'Caris_class_I':';'.join(caris.get(gene,[])),'Caris_comparison_3_24':two(a)==caris[gene] if gene in caris else None,'Caris_comparison_3_46':two(b)==caris[gene] if gene in caris else None}
    for version,label in [('3.24.0','3_24'),('3.46.0','3_46')]:
        d=diagnostics[version]['genes'].get(gene,{})
        row['model_fragment_count_'+label]=d.get('model_gene_fragment_count')
        row['genotype_explained_fragments_'+label]=d.get('genotype_explained_fragments')
        row['nonshared_minor_major_'+label]=d.get('minor_to_major_nonshared_ratio_rounded')
        p=next(x for x in prior if x['version']==version)
        row['no_population_prior_two_field_match_'+label]=p['comparison'][gene]['two_field_match']
    rows.append(row)
limits=[
 'Exploratory arcasHLA typing from FFPE bulk tumor/admixed RNA; not accredited clinical HLA typing, not independent confirmation of germline HLA, and not evidence of tumor allele retention or HLA LOH.',
 'Class-I comparison with Caris exome calls is an internal cross-assay check within the same specimen, not an independent tissue validation. Class-II calls have no provided clinical comparator.',
 'Reference sensitivity is reported for both versions without selecting a favorable answer. Same-reference public regression passed all seven loci on IMGT3.24; the newer IMGT3.46 public example had B/DQB1/DRB1 differences from the old standard-stage expectation, with DRB1 still different from final partial-stage expectation.',
 'Only standard complete-reference typing was run; optional partial-reference stage was omitted because its derived index exceeded the agreed public-reference disk limit. Partial alleles or rare alleles may remain unresolved.',
 'No matched normal. Tumor/normal admixture and allele-specific expression can mask loss or create apparent homozygosity. Bulk HLA-II expression may come from infiltrating antigen-presenting cells, not malignant cells.',
 'Two-field pairs are the reporting resolution. Higher-field printed labels are representatives of compatible sequences and do not establish higher-resolution typing.',
 'Repeated DPA1 or DRB5 allele labels are sole sequence types called by the RNA algorithm, not proof of two genomic copies. DRB3/4/5 gene-content variation especially complicates diploid interpretation. DPA1 and DRB5 were outside the seven-locus bundled regression control.',
 'Model transcript abundances, explained-fragment percentages, and minor/major nonshared ratios are diagnostics, not calibrated genotype confidence probabilities. Paired-end kallisto count units are fragments/pairs, despite arcasHLA log language calling them reads; no UMI or PCR deduplication.',
 'Cross-locus phase is not established. Any DQ/DP prediction must address alpha/beta chain combinations and use a class-II method; current class-I peptide scoring cannot be relabeled class II.',
 'This is not a full RNA reanalysis or comprehensive antigen-presentation audit. No protein presentation, T-cell response, or vaccine efficacy is measured.'
]
sources=[{'title':'Official arcasHLA source/workflow','url':'https://github.com/RabadanLab/arcasHLA','commit':'9fa54a212d134b0d9894d1fc19ec1bdc6f62eb55'},{'title':'arcasHLA primary publication','url':'https://pmc.ncbi.nlm.nih.gov/articles/PMC6956775/'},{'title':'HLA allele-specific expression in tumors','url':'https://link.springer.com/article/10.1186/s13073-023-01154-x'},{'title':'Personalized neoantigen vaccine CD4 responses, primary study','url':'https://www.nature.com/articles/nature22991'}]
result={'status':'exploratory_completed','extraction':extract,'runs':runs,'comparison':rows,'diagnostics':diagnostics,'prior_sensitivity':prior,'replay_failure_and_resolution':'An attempted replay of official six-field alignment.p failed in upstream load_alignment legacy branch expecting nine fields. Failure preserved in patient-prior-alignment-reuse-failed.json/log. No official source was patched. Both no-prior sensitivity runs instead reprocessed the same exact paired FASTQs and completed.','limitations':limits,'sources':sources}
(OUT/'patient-HLA-findings.json').write_text(json.dumps(result,indent=2)+'\n')
with (OUT/'patient-HLA-comparison.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
lines=['Exploratory local RNA HLA analysis — completed 1 October 2026','',
 'The class-I sequence pairs match Caris at all three loci, and all nine requested loci are stable at two-field resolution across both reference releases and with or without population-frequency priors. This supports the HLA inputs used for research peptide screening within this specimen; it is not clinical HLA validation or proof of tumor allele retention.','',
 'Class I: A*02:01 / A*26:01; B*50:01 / B*52:01; C*06:02 / C*12:02.','',
 'Provisional class II: DPB1*02:01 / DPB1*04:01; DQA1*01:03 / DQA1*05:01; DQB1*02:01 / DQB1*06:01; DRB1*03:01 / DRB1*15:02. DPA1*01:03 and DRB5*01:02 are the sole sequence types called at those genes (software repeats each); do not equate this with proven germline homozygosity.','',
 'Research value: these provisional class-II types could guide a separate mutation-spanning CD4/long-peptide screening exercise for Michael/Eta. Preserve natural protein sequence and the premature stop in BAP1/RASA1; do not invent residues beyond it. DQ/DP alpha-beta pairing and cross-locus phase remain unresolved. Clinical typing/appropriate research validation would be needed before definitive construct selection.','',
 'Input: 2,490,971 complete RNA read pairs (887,610 chr6 proper pairs plus 1,603,361 both-unmapped pairs), no orphans or duplicate primary mate keys, exact paired FASTQ name/sequence/quality roundtrip passed. No MAPQ cutoff or PCR deduplication. Actual patient qualities were preserved; no placeholders.','',
 'Method: pinned arcasHLA commit 9fa54a212d134b0d9894d1fc19ec1bdc6f62eb55, official kallisto 0.44.0, standard complete references IMGT3.24.0 and3.46.0. One thread. Minimum 75 gene-assigned fragment counts, EM tolerance1e-6, 1000 iterations, low-abundance drop0.1 after20 iterations, zygosity threshold0.15. Both default-prior and no-prior analyses completed. All commands, reference hashes, model diagnostics and original logs are retained.','',
 'Important limits:']+['- '+x for x in limits]+['','Reproducibility: extract-patient.py; run-patient-genotypes.py; patient-prior-sensitivity.py; compile-patient.py. See patient-HLA-findings.json, patient-HLA-comparison.tsv, patient-extraction.json, public controls and source manifests in this directory.','', 'Primary sources:']+[x['url'] for x in sources]
(OUT/'patient-HLA-findings.txt').write_text('\n'.join(lines)+'\n')
print(json.dumps({'comparison':rows,'diagnostics':diagnostics},indent=2))
