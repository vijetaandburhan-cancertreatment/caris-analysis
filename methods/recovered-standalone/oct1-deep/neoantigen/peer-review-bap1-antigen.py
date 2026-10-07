"""Independent bounded audit: reference translations, exhaustive neighbors, tables.

No patient sequences uploaded. Recomputes normal neighbors using every window
and vectorized Hamming distances, independent of the production seed search.
"""
from pathlib import Path
from collections import Counter
import csv, gzip, hashlib, json, time
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

OUT = Path(__file__).resolve().parent
started = time.time()
context = json.loads((OUT/'BAP1-RASA1-local-protein-context.json').read_text())
translations = {}
for gene, accession, start, stop in [('BAP1','NM_004656.3',168,178),('RASA1','NM_002890.2',747,747)]:
    r = SeqIO.read(OUT/'reference-transcripts'/f'{accession}.gb', 'genbank')
    cds = next(f for f in r.features if f.type == 'CDS')
    sequence = str(cds.extract(r.seq))
    wt = str(Seq(sequence).translate(cds=True))
    assert wt == cds.qualifiers['translation'][0]
    edited = sequence[:start-1] + sequence[stop:]
    mutant = str(Seq(edited[:len(edited)//3*3]).translate(to_stop=True))
    first = next(i for i,(a,b) in enumerate(zip(wt,mutant)) if a != b)
    matches = []
    for c in context['records']:
        if c['gene'] != gene: continue
        protein = wt if c['sequence_type']=='wt' else mutant
        seq = protein[c['amino_acid_start1']-1:c['amino_acid_end1']]
        assert seq == c['sequence']
        assert hashlib.sha256(seq.encode()).hexdigest() == c['sha256_sequence']
        assert c['first_altered_position1'] == first+1
        matches.append(c['sequence_type'])
    translations[gene] = {'accession': accession, 'deleted_coding_bases': sequence[start-1:stop], 'first_changed_aa1': first+1, 'mutant_protein_length': len(mutant), 'altered_tail': mutant[first:], 'context_types_reproduced':matches}
    if gene == 'BAP1':
        translations[gene]['peptides'] = []
        for peptide in ['IEERKGLYL','EERKGLYL']:
            position = mutant.index(peptide)
            translations[gene]['peptides'].append({'peptide':peptide,'start_aa1':position+1,'end_aa1':position+len(peptide),'same_position_wildtype':wt[position:position+len(peptide)],'remaining_residues_before_stop':mutant[position+len(peptide):]})

queries = ['IEERKGLYL','EERKGLYL']
recomputed = Counter()
records = residues = 0
with gzip.open(OUT/'gencode.v37.pc_translations.fa.gz','rt') as handle:
    for rec in SeqIO.parse(handle,'fasta'):
        sequence = str(rec.seq)
        residues += len(sequence); records += 1
        values = np.frombuffer(sequence.encode(),dtype=np.uint8)
        for query in queries:
            n = len(query)
            if len(sequence)<n: continue
            # Every legal start, no seeds and no early candidate restriction.
            distances = np.zeros(len(sequence)-n+1,dtype=np.uint8)
            for i, letter in enumerate(query.encode()):
                distances += values[i:i+len(distances)] != letter
            for offset in np.flatnonzero(distances<=2):
                offset = int(offset)
                recomputed[(query,sequence[offset:offset+n],rec.description,offset+1)] += 1
expected = Counter()
neighbors = json.loads((OUT/'bap1-normal-sequence-neighbors.json').read_text())
for row in neighbors['results']:
    for occurrence in row['occurrences']:
        expected[(row['query'],row['normal_reference_peptide'],occurrence['reference_header'],occurrence['start_aa1'])] += 1
assert recomputed == expected, {'missing':len(expected-recomputed),'extra':len(recomputed-expected)}
neighbor_summary = {}
for q in queries:
    unique = {key[1] for key in recomputed if key[0]==q}
    neighbor_summary[q] = {'unique_by_substitution_count':{str(d):sum(sum(a!=b for a,b in zip(q,s))==d for s in unique) for d in range(3)},'all_occurrences':sum(v for k,v in recomputed.items() if k[0]==q)}

affinity = list(csv.DictReader((OUT/'bap1-normal-neighbor-affinity.tsv').open(),delimiter='\t'))
assert len(affinity)==56
assert all(r['allele']=='HLA-B*50:01' for r in affinity)
checks=[]
for release, filename in [('2.3.0','mhcflurry-focused-affinity-validation.csv'),('2.2.0','mhcflurry-historical-affinity-check.csv')]:
    reference = list(csv.DictReader((OUT/filename).open()))
    for row in reference:
        if row['allele']!='HLA-B*50:01': continue
        match = next((r for r in affinity if r['model_release']==release and r['peptide']==row['peptide']),None)
        if not match: continue
        value = float(match['predicted_affinity_nM']); expected_value = float(row['mhcflurry_affinity'])
        assert abs(value-expected_value) <= max(0.01,1e-5*expected_value)
        checks.append({'release':release,'peptide':row['peptide'],'relative_difference':abs(value-expected_value)/expected_value})

result = {'status':'PASS with presentation clarifications','elapsed_seconds':time.time()-started,'sequence_context_reconstruction':translations,'normal_neighbor_recount':{'implementation':'Unseeded vectorized Hamming distance for every equal-length window in every reference record','records':records,'amino_acid_residues':residues,'exact_full_occurrence_multiset_matches':True,'summary':neighbor_summary},'fixed_allele_affinity_summary':{'rows':len(affinity),'allele':'HLA-B*50:01','cross_table_checks':checks,'new_model_inference_run':False},'clarifications':[{'artifact':'candidate-screen-summary.tsv','issue':'RNA_alt_fragments is variant-level support, not full peptide-window support. BAP1 254 describes the deletion; focused full-window RNA support is 241 for IEERKGLYL and 245 for EERKGLYL.','suggestion':'Rename column RNA_variant_alt_fragments or add explicit legend; use focused counts when claiming complete local peptide sequence support.'},{'artifact':'candidate-screen-summary.tsv','issue':'Mutant and wild-type best alleles are selected independently. Their displayed best-affinity values cannot always be used as same-HLA mutation-effect comparisons. The focused BAP1 B*50:01 comparison is correctly allele matched.','suggestion':'Keep allele labels and explicitly prohibit cross-allele affinity-ratio interpretations.'},{'artifact':'bap1-antigen-interpretation.txt','issue':'Current tumor confirmation is a research development decision and need not automatically mean a new biopsy.','suggestion':'Prefer existing suitable tumor material where available; any new biopsy remains a clinical decision.'}],'interpretation_limits':['Same-data independent calculations are computational reproducibility checks, not independent clinical validation.','Two overlapping peptides are one mutation hypothesis; two MHCflurry model releases are correlated model checks.','The best predictions were selected after a broader 676-peptide screen and carry no calibrated patient-specific success probability or false-discovery rate.','Predicted affinity/presentation scores do not establish binding measurements, malignant-cell presentation, immune recognition, response, safety or survival.','No matched normal establishes somatic status; GENCODE normal-reference absence and Hamming distance do not establish tumor specificity or TCR cross-reactivity.','RNA supports the local sequence; it does not establish full-length mutant isoform, stable protein, NMD escape rate or HLA retention.']}
(OUT/'peer-review-bap1-antigen.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'elapsed_seconds':result['elapsed_seconds'],'neighbor_recount':result['normal_neighbor_recount'],'affinity_cross_table_checks':len(checks)},indent=2))
