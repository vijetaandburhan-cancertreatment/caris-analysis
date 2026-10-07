from pathlib import Path
import csv,json,hashlib,time
OUT=Path(__file__).resolve().parent
x=json.loads((OUT/'results.json').read_text());r={v['kind']:v for v in x['results']}
rows=[]
for kind,v in r.items():
    for z in v['peptide_haplotype_support']:
        a='exact_deletion_haplotype'
        rows.append({'kind':kind,'peptide_or_context':z['label'],'exact_mutant_nucleotide_reads':z['reads'].get(a,0),'exact_mutant_nucleotide_query_name_fragments':z['fragments'].get(a,0),'exact_reference_query_name_fragments':z['fragments'].get('exact_reference_haplotype',0),'other_haplotype_fragments':z['fragments'].get('other_haplotype',0),'discordant_fragments':z['fragments'].get('discordant',0),'exact_mutant_forward_reads':z['read_strands'].get(a+'_forward',0),'exact_mutant_reverse_reads':z['read_strands'].get(a+'_reverse',0),'distinct_mutant_alignment_patterns':z['distinct_alignment_patterns'].get(a,0),'clean_mutant_fragments':z['clean_no_softclip_at_least_5bp_from_ends_fragments'].get(a,0),'interpretation':'Locally phased nucleotides; conditional translation, not measured peptide/protein; fragments not independent UMI molecules.'})
with (OUT/'BAP1-peptide-local-evidence.tsv').open('w') as f:
    w=csv.DictWriter(f,rows[0].keys(),delimiter='\t');w.writeheader();w.writerows(rows)
findings={
    'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
    'conclusion':'The exact local nucleotide haplotypes that conditionally translate to IEERKGLYL and EERKGLYL are strongly supported by both DNA and RNA. No credible nearby linked alteration found in this bounded review contradicts those dominant local sequences. This is nucleotide/short-read validation, not evidence of protein production or HLA presentation.',
    'exact_haplotype_evidence':rows,
    'reference_and_mapping':{'genomic_event':'hg38 chr3:52408550 CGCCGGGACCGG>C','RefSeq':'NM_004656.3','coding_event':'c.168_178del11','conditional_protein':'p.(Arg57LysfsTer8)','MANE_GENCODE37':'ENST00000460680.6','affected_exon':4,'affected_exon_genomic_0based_halfopen':[52408473,52408606],'affected_exon_reference_coding_positions_1based':[123,255],'9mer_reference_coding_window':[157,194],'9mer_mutant_coding_window':[157,183],'9mer_genomic_window0':[52408534,52408572],'8mer_reference_coding_window':[160,194],'8mer_genomic_window0':[52408534,52408569],'extended_context':'FKWIEERKGLYLGG*','extended_reference_coding_window':[148,203],'target_exon_public_hg38_RefSeq_exact_match':True},
    'nearby_variants':{
        'scope':'375 coding bases in MANE exons 1--5; all qualifying base observations retained in nearby-coding-base-evidence.tsv.',
        'supplied_VCF':'Only the reported 11bp BAP1 deletion is present in the bounded five-exon region.',
        'DNA':'No additional SNV reaches >=5 non-discordant query-name fragments and >=1%. Target-exon strict DNA read depth spans 1216--2696x. A 1bp deletion at c.210 is seen in six query-name fragments but every supporting read is within 5bp of its alignment end; it is after the predicted early stop and is not promoted.',
        'RNA_peptide_window':'No SNV in the exact nine-/eight-peptide coding window reaches >=3 non-discordant query-name fragments and >=1%. The dominant complete haplotypes match the predicted sequence directly.',
        'RNA_upstream_flank':'At c.154 (immediately upstream of the nine-peptide window), 17 G and 8 C forward-genomic alternate fragments are unphased to the known deletion; all fail the no-softclip/at-least-5bp-from-end check. They do not contradict the 231 exact full-context RNA fragments spanning c.148--203.',
        'RNA_downstream':'Four other low-level RNA SNV observations are at reference c.265,289,312,337, all beyond the modeled stop and outside the candidate peptide window. Some have poor end/clipping or repeated-coordinate support; these are retained as minor observations, not validated variants.',
        'RNA_minor_indels':'A 12bp deletion overlapping the known 11bp event has five reads/four query-name fragments with exact local sequence; no exact 12bp haplotype in the corresponding stringent DNA check. This rare RNA observation is unresolved and not silently discarded. The dominant 11bp RNA haplotype has 410 reads/241 fragments over the same nine-peptide window. Other small RNA indel observations lack phase to the main deletion and do not change the locally established sequence.',
        'meaning':'A clinical negative or absence of every subclonal event is not established. The local peptide sequence can coexist with minor other RNA molecules or technical errors.'
    },
    'splice_evidence':{
        'exon3_to4':next(z for z in r['RNA']['junctions'] if z['annotation']=='exon_3_to_4'),
        'exon4_to5':next(z for z in r['RNA']['junctions'] if z['annotation']=='exon_4_to_5'),
        'both_junctions':r['RNA']['dual_exon4_junction_support'],
        'interpretation':'Normal exon-4 inclusion is directly supported on deletion-bearing local RNA fragments. Twenty-four paired query-name fragments support the deletion and both neighboring junctions. This does not identify a full-length isoform or establish reference translational initiation.',
        'minor_junctions':'Two fragments support the exon3-to5 junction skipping exon4; one fragment each supports two other local junctions. Skipping fragments cannot be phased to a deletion in the omitted exon. Counts are not an isoform abundance/PSI estimate.'
    },
    'research_value':'Keep these BAP1 peptides on the conditional vaccine research shortlist with strengthened local sequence provenance. Do not label them verified antigens or clinically suitable vaccine targets; tumor specificity, actual protein/peptide production, HLA presentation and recognition remain separate questions.',
    'limitations':x['limitations'],
    'files':{'primary_script':'validate.py','primary_result':'results.json','reference':'reference-provenance.json','local_haplotype_table':'BAP1-peptide-local-evidence.tsv','minor_observation_audit':'minor-observation-audit.json','minor_12bp_check':'minor-12bp-deletion-haplotype-check.json','splicing':'exon4-splice-junction-evidence.tsv'},
    'code_SHA256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(OUT.glob('*.py'))}
}
(OUT/'findings.json').write_text(json.dumps(findings,indent=2)+'\n')
text='''The BAP1 local sequence check supports both proposed peptide sequences.

Exact mutant nucleotide haplotypes spanning the full candidate window:
                          DNA fragments   RNA fragments
IEERKGLYL                       691              241
EERKGLYL                        693              245
FKWIEERKGLYLGG* context          594              231

Every counted supporting read spans every nucleotide in its window and the
deletion at Q20 or better. Mates are collapsed by read group + query name.
These are not independent UMI molecules. RNA exact nine-peptide support has
206 forward / 204 reverse reads and 218 distinct alignment-coordinate patterns;
192 fragments remain without soft clipping and at least 5bp from read ends.

The expected exon3-to4 junction has 484 fragments (143 phased to the deletion),
and exon4-to5 has 249 (52 phased). Twenty-four deletion-bearing paired fragments
jointly support both neighboring junctions. This supports local exon inclusion,
not a full-length translated isoform. Rare exon4-skipping support (two fragments)
is also recorded; these counts are not a quantitative isoform fraction.

No additional DNA SNV meets the prespecified >=5 fragments / >=1% review
threshold over coding exons1--5. No RNA SNV meets >=3 / >=1% within the actual
peptide window. Minor RNA mismatches immediately upstream fail stringent
end/clipping checks; more distant mismatches are outside the peptide/after the
modeled stop. Four RNA fragments have a 12bp rather than 11bp deletion, versus
241 with the exact 11bp event in the same window; the minor 12bp haplotype has
no exact DNA support in that stringent check. It remains unresolved, not a
validated separate allele or a reason to replace the dominant peptide sequence.

Thus nearby linked sequence does not contradict the dominant proposed local
translation. The defensible wording is 'read-supported local nucleotide
sequence, conditionally encoding these peptides.' We have not measured peptide
production, HLA presentation, T-cell recognition, tumor specificity, full-length
transcript phase, translational initiation or an NMD probability.
'''
(OUT/'findings.txt').write_text(text)
print(text)
