"""Summarize existing MHCflurry outputs without rerunning predictions.

Preserves numeric output strings. Mutant/WT best alleles are separately chosen
by the six-allele model call and cannot generally be used for matched-HLA ratios.
"""
from pathlib import Path
import csv, json, hashlib

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
target = ROOT/'outputs/caris-deep-review'
target.mkdir(parents=True,exist_ok=True)
peptides = {r['candidate_id']:r for r in csv.DictReader((OUT/'candidate-peptides.tsv').open(),delimiter='\t')}
models = list(csv.DictReader((OUT/'mhcflurry-predictions.csv').open()))
combined = {(r['candidate_id'],r['sequence_type']):r for r in models if ';' in r['allele']}
assert len(combined)==2*len(peptides)
rows=[]
for candidate,source in peptides.items():
    mt = combined[(candidate,'mutant')];wt = combined[(candidate,'wildtype')]
    row={k:mt[k] for k in ['candidate_id','variant','gene','peptide','best_allele','affinity','affinity_percentile','processing_score','presentation_score','presentation_percentile']}
    row.update({'wildtype_peptide':wt['peptide'],'wildtype_best_allele':wt['best_allele'],'wildtype_best_affinity_nM':wt['affinity'],'normal_reference_exact_match_count':source['normal_reference_exact_match_count'],'RNA_variant_alt_fragments':source['RNA_alt_fragments'],'clinical_status':'unvalidated candidate; somatic status unknown'})
    rows.append(row)
rows.sort(key=lambda r:float(r['presentation_percentile']))

# Check the first revised generation against the pre-existing table: no model
# value, sequence, HLA choice or support count may change under this label edit.
old_path=OUT/'candidate-screen-summary.tsv'
comparison = {'preexisting_rows':None,'all_data_preserved':None}
if old_path.exists():
    old=list(csv.DictReader(old_path.open(),delimiter='\t'))
    for row in old:
        if 'RNA_alt_fragments' in row:
            row['RNA_variant_alt_fragments']=row.pop('RNA_alt_fragments')
    assert len(rows)==len(old)
    assert rows==old, 'Existing summary data or ordering differs from reproducible source reconstruction'
    comparison={'preexisting_rows':len(old),'all_data_preserved':True}

for path in [old_path,target/'Peptide-research-screen.tsv']:
    with path.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(rows)

schema='''Peptide research screen — interpretation and column definitions

This is a research model screen, not a list of validated neoantigens or a vaccine
design. It preserves all 676 candidate rows from the selected 21 source variants
(19 produce altered-amino-acid windows in this translation model). Somatic
status is unestablished without a matched normal comparator. Several alleles
occur in population catalogs. Exact absence from the selected GENCODE37 normal
protein reference does not establish tumor specificity or safety.

RNA_variant_alt_fragments counts paired-read identifiers supporting the variant
locus, after the documented allele-evidence filters. It does not count complete
peptide-encoding windows, independent UMI molecules, or translated proteins.
For BAP1, 254 is the variant-level RNA count; the separate focused local audit
found 241 RNA paired names supporting the complete IEERKGLYL nucleotide window
and 245 supporting the complete EERKGLYL window. Do not interchange those counts.

best_allele/affinity are the mutant peptide's minimum predicted affinity across
the six supplied HLA-I alleles. wildtype_best_allele/wildtype_best_affinity_nM
are optimized independently for its corresponding wild-type peptide. A ratio
between these displayed affinities is NOT generally a matched-HLA effect.
Use the full per-allele mhcflurry-predictions.csv for a fixed-allele comparison.
The focused BAP1 HLA-B*50:01 mutant/WT comparison is allele matched.

affinity is predicted binding IC50 in nM, not a laboratory measurement.
affinity_percentile and presentation_percentile are model reference-distribution
rankings, not chances of presentation, response or clinical benefit.
processing_score and presentation_score are model outputs, not measurements.
The row sort uses combined-six-allele presentation_percentile. Selected best
scores after a broad screen do not have a calibrated false-discovery rate or
patient-specific success probability. Two MHCflurry releases are correlated
model checks; two overlapping BAP1 peptides represent one mutation hypothesis.

normal_reference_exact_match_count counts GENCODE37 protein records containing
the peptide, not tissue expression or every within-protein occurrence. That
reference omits normal population polymorphisms and noncanonical ORFs.

Source generation: generate-candidate-summary.py from candidate-peptides.tsv
and mhcflurry-predictions.csv; no predictions are rerun during summarization.
Sequence assumptions, versions and model hashes are in candidate-preparation.json
and mhcflurry-run-manifest.json. Clinical benefit, antigen presentation, immune
recognition, HLA retention and normal-cell safety require separate evidence.
'''
(OUT/'candidate-screen-summary-schema.txt').write_text(schema)
(target/'Peptide-research-screen-schema.txt').write_text(schema)
manifest={'generator':str(Path(__file__).relative_to(ROOT)),'candidate_rows':len(rows),'verification':comparison,'RNA_column_renamed_from':'RNA_alt_fragments','RNA_column_renamed_to':'RNA_variant_alt_fragments','source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [OUT/'candidate-peptides.tsv',OUT/'mhcflurry-predictions.csv']},'output_sha256':hashlib.sha256(old_path.read_bytes()).hexdigest()}
(OUT/'candidate-screen-summary-generation.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({'candidate_rows':len(rows),'verification':comparison,'outputs_identical':old_path.read_bytes()==(target/'Peptide-research-screen.tsv').read_bytes()},indent=2))
