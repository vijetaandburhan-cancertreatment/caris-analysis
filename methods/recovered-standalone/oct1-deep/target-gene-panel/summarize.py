"""Summarize the completed bounded panel without re-running BAM analysis."""
from pathlib import Path
import json,csv,hashlib,collections
P=Path(__file__).resolve().parent
j=json.loads((P/'screen-results.json').read_text()); setup=json.loads((P/'setup.json').read_text())
cs=j['candidates']; gs=j['gene_summaries']
assert len(cs)==54 and all(c['Caris_VCF_exact_matches'] for c in cs)
assert sum(x['CDS_plus4_bases'] for x in gs)==87057
assert all(x['CDS_plus4_bases']==x['bases_ge20'] for x in gs)
protein=[c for c in cs if c['consequence'] in ('missense','inframe_indel','frameshift','stop_gain')]
assert len(protein)==18
labels=collections.Counter(m['clinical'] for c in protein for m in c['Caris_VCF_exact_matches'])
assert labels=={'Benign':17,'VUS':1}
findings={
 'status':'Completed bounded research screen; no clinical variant reclassification',
 'main_finding':'All 54 candidates already occur as exact alleles in the supplied Caris VCF export. This 22-gene selected-CDS screen found no additional protein-altering candidate above its filters.',
 'genes':setup['genes'],
 'reference_and_transcripts':'hg38 generic public reference windows; selected GENCODE37 coding transcripts: 21 MANE Select plus NTRK3 ENST00000394480.6 APPRIS principal_1. Whole transcript isoform repertoire and regulatory regions not assessed.',
 'coverage':{'selected_CDS_plus4_positions':87057,'positions_with_at_least20_strict_quality_reads':87057,'minimum_over_all_selected_positions':min(x['minimum_strict_read_depth'] for x in gs),'interpretation':'Observed qualifying base-read depth, not validated diagnostic sensitivity or a comprehensive negative result.'},
 'counts':dict(collections.Counter(c['consequence'] for c in cs)),
 'exact_source_matches':54,
 'protein_changing_source_labels':dict(labels),
 'source_label_limit':'Benign and VUS here are Caris labels, not independent clinical reclassifications. Variant fraction does not establish somatic or germline status without a matched normal.',
 'existing_ALK_control':{'hg38':'chr2:29222407 G>A','selected_transcript':'ENST00000389048.8','conditional_protein':'T1151M','Caris_label':'VUS','initial_strict_DNA_alt_ref_paired_names':[791,1436],'initial_strict_RNA_alt_ref_paired_names':[41,102],'count_source':'outputs/caris-analysis/Caris-candidate-evidence.tsv','audit_limit':'Existing primary-script strict recount; the table does not claim a second independent recount for ALK. Not a fusion, clinical driver or treatment-sensitivity call.'},
 'FGFR1_representation':{'hg38':'chr8:38428395 G>GTCA','screen_consequence':'in-frame insertion; precise HGVS not assigned by screen','Caris_protein_label':'D133dup','Caris_classification':'Benign','note':'One added Asp in an existing poly-Asp tract. Original S134D first-difference label was misleading as HGVS and removed; exact genomic allele/counts retained. Repeat-shifted CIGAR placement means narrow-window RNA haplotype counts cannot be used to infer absent expression. Independent peer audit provides longer-haplotype context.'},
 'ROS1_numbering':'Selected transcript protein residue numbers differ from Caris transcript numbering for five missense records. The exact hg38 allele matches are preserved; this is not evidence that their genomic alleles disagree.',
 'POLE_flank_events':[{'hg38':'chr12:132661167 G>GA','Caris_line':2164,'Caris_filter':'sb;R8.1'},{'hg38':'chr12:132661167 GA>G','Caris_line':2166,'Caris_filter':'R8.1'}],
 'POLE_interpretation':'Both are previously filtered source alleles near a coding-exon boundary, not coding protein alterations in this selected transcript. They remain unpromoted. Repeat alignment/indel representation and possible splice context require more than this screen; RNA-window nonassessability is not a negative splice finding.',
 'method':j['method'],
 'limits':[
  'Bounded candidate generation, not a full exome SNV re-call, comprehensive structural-variant or splice analysis, or clinical testing.',
  'No matched normal, calibrated purity/copy-number model, or verified assay capture BED; no somatic, biallelic-loss, germline or clinical negative claims.',
  'Single selected coding transcript per gene and single-variant translation. Stop codons outside the GENCODE CDS, other isoforms, deep intronic/regulatory alterations and complex phased alleles are not comprehensively evaluated.',
  'Depth is qualified read depth; PCR/template independence is not established. Query-name collapse for indel screening is not UMI validation.',
  'Screen SNVs require at least 10 alternate reads and 5% allele fraction; initial indels require at least 10 query-name fragments and use alternate-read/anchor-read depth. Counts at the screening stage are not harmonized clinical VAFs; indels were not normalized.',
  'Good coding coverage does not exclude kinase fusions, copy-number changes or clinically important variants below the chosen thresholds.'
 ],
 'runtime_seconds':j['elapsed_seconds'],
 'resource_note':'Guard observed a sampled peak RSS of 55,869,440 bytes over the last ~1 second; this is not a measured full-run peak. Full screen completed in 52.37 seconds; no heavy or ongoing job remains.',
 'independent_audit':'Origin peer independently reconstructed all 22 reference proteins and 52 coding consequences against GENCODE37 protein FASTA. Additional exact public-allele/read-context audit is saved separately; see audit paths supplied in final handoff.',
 'source_sha256':{x:hashlib.sha256((P/x).read_bytes()).hexdigest() for x in ['screen.py','screen-results.json','target-annotation.json','setup.json']},
}
audit_paths=['audit.txt','audit.json','independent-audit.py','independent-audit-data.json','independent-public-audit.py','independent-public-audit.json']
assert all((P/x).exists() for x in audit_paths)
findings['independent_audit']='PASS after the incorporated transcript/indel label corrections; all 22 reference proteins, 52 coding edits and 54 exact VCF matches independently checked. Full repeat-sensitive FGFR1/POLE recounts and public allele context are in the audit files.'
findings['independent_audit_sha256']={x:hashlib.sha256((P/x).read_bytes()).hexdigest() for x in audit_paths}
findings['public_allele_context']={'exact_missense_catalog_matches':17,'missense_with_global_reference_frequency_at_least_1percent':16,'ALK_T1151M':{'id':'rs113994091','gnomADe_alt_frequency':0.00001642,'gnomADg_alt_frequency':0.00004602},'limit':'Generic public allele queries only, no patient/raw sequence upload. Public variation-level clinical metadata can mix alleles/conditions and was not treated as allele-specific clinical adjudication. Population occurrence does not identify constitutional status for this patient or establish a treatment target.'}
(P/'findings.json').write_text(json.dumps(findings,indent=2)+'\n')
with (P/'coverage.tsv').open('w') as f:
 w=csv.DictWriter(f,list(gs[0]),delimiter='\t');w.writeheader();w.writerows(gs)
fields=['gene','chrom','pos1','ref','alt','consequence','protein','transcript','Caris_protein','Caris_classification','Caris_filter','Caris_VCF_line']
with (P/'protein-changing-source-context.tsv').open('w') as f:
 w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader()
 for c in protein:
  m=c['Caris_VCF_exact_matches'][0]
  w.writerow({**{k:c.get(k,'') for k in fields[:8]},'Caris_protein':m['protein'],'Caris_classification':m['clinical'],'Caris_filter':m['filter'],'Caris_VCF_line':m['line']})
text='''22-gene coding-screen findings\n\nAll 54 candidates already occur as exact alleles in the supplied Caris VCF. This bounded screen found no additional protein-altering candidate above its filters.\n\n87,057 selected CDS plus/minus 4-bp positions had at least 20 qualifying reads (minimum 146). Transcripts were 21 GENCODE37 MANE Select plus NTRK3 APPRIS principal_1. This is observed depth, not a clinically validated negative result.\n\nThe 54 candidates comprise 34 synonymous, 17 missense, one in-frame insertion and two noncoding/CDS-flank events. Of the 18 protein-changing candidates, Caris labels 17 Benign and the existing ALK T1151M VUS. These are source labels, not independent clinical reclassifications. No new treatment target is established.\n\nFGFR1 is an in-frame insertion in a poly-Asp repeat; the screen does not assign precise HGVS. Its old first-difference S134D label was removed because it could be misread as a missense. Caris reports D133dup, Benign. RNA CIGAR placements shift in the repeat; only longer exact haplotypes are interpretable.\n\nBoth POLE exon-flank insertion/deletion events were already filtered in the Caris VCF (sb;R8.1 and R8.1) and are not promoted. RNA-window nonassessability cannot establish absence of an expressed variant or normal splicing.\n\nThe companion JSON supplies filters, limitations, hashes and existing ALK read controls. This is not a full-exome SNV re-call, a comprehensive splice/structural exclusion, proof of tumor-specificity, or clinical validation. Per-gene depth is in coverage.tsv; genomic alleles and source labels are in protein-changing-source-context.tsv.\n'''
text += '\nIndependent audit completed: audit.txt/audit.json document matching reference proteins, coding edits, exact source alleles and repeat-sensitive read context. All 17 missense alleles have exact public catalog mappings; 16 occur at >=1% in global reference data. ALK T1151M is rare in public data and remains a source VUS, not a proven somatic or actionable call. No additional material correction remains.\n'
(P/'findings.txt').write_text(text)
print(json.dumps({'candidates':len(cs),'protein_changing':len(protein),'classifications':dict(labels),'positions':87057,'outputs':['findings.json','findings.txt','coverage.tsv','protein-changing-source-context.tsv']}))
