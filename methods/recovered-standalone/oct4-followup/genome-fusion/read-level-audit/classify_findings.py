#!/usr/bin/env python3
"""Manual scientific classification layered on reproducible exact-match evidence; no clinical calls."""
import collections,csv,hashlib,json,pathlib
base=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct4-fusion-read-audit')
p=json.load(open(base/'patient-evidence/candidate-evidence.json'));s={r['candidate_id']:r for r in json.load(open(base/'supplemental-evidence/candidate-evidence.json'))};pa=json.load(open(base/'patient/junction-read-audit.json'));sa=json.load(open(base/'supplemental/junction-read-audit.json'));meta={r['candidate_id']:r for r in pa['rows']};sm={r['candidate_id']:r for r in sa['rows']};shared=json.load(open(base/'patient-evidence/shared-marker-groups.json'))
manual={
0:('credible_DNA_RNA_sequence_event_mechanism_unresolved','Abundant diverse exact RNA junction support plus100 DNA fragments with an exact100nt marker (55Q30), independently recounted. Self-event joins 3prime UTR to intronic sequence, with no caller peptide; not an intergenic driver fusion or demonstrated coding neoantigen. RNA processing/circularization, DNA change, normal-cell contribution and germline status remain unresolved.'),
1:('credible_DNA_RNA_repeat_insertion_somatic_status_unknown','Abundant RNA and DNA support for a repeat-rich27bp insertion representation; exact100nt DNA marker329Q20/262Q30 fragments independently recounted. Caller duplication plus local mismatch has an equivalent insertion representation. Repeat-alignment mate conflicts preclude a reliable VAF; neither tumor specificity nor driver/neoantigen function is established.'),
2:('artifact_compatible_repeat_nonunique_marker','GC repeat junction markers occur in many ordinary genome/transcript positions. Counts are not specific evidence for a NOTCH2NLC ITD; gene assignment/repeat length remains unresolved.'),
3:('repeat_or_paralog_confounded','One-base insertion consensus is raw-read supported, but both anchors have multiple genomic matches and most caller names have alternative STAR chimeric alignments. Same canonical marker as the other chr15 row; not two independently established inversions.'),
4:('repeat_or_paralog_confounded','Reverse-complement-equivalent one-base insertion marker to the other chr15 row. Multiple loci per anchor and alternative STAR alignments prevent a unique rearrangement interpretation.'),
5:('repeat_or_paralog_confounded','Exact raw sequence exists, but one anchor has thousands of genomic matches and all20 caller-named pairs have alternative STAR chimeric alignments. 3prime-3prime geometry does not establish a conventional chimeric protein.'),
6:('repeat_or_paralog_confounded','Identical exact marker is assigned to three ACE partner calls. Partner anchor has extreme genome repetition; sequence corroboration does not identify a specific partner or DNA rearrangement.'),
7:('RNA_supported_nonproductive_self_event_unresolved','Exact raw support for an intragenic/intronic inverted3prime-3prime junction. No conventional fusion peptide; library hairpin/RNA-processing versus genomic inversion remains unresolved.'),
8:('low_support_nonproductive_self_event_unresolved','Three strict50nt/Q20 pairs for an out-of-frame inverted self-junction; oneQ30 pair. Too little evidence to establish a biological protein product or tumor-specific inversion.'),
9:('repeat_or_paralog_confounded','Insertion consensus is raw-supported but PRAME-family/repetitive partner anchors are nonunique, and all caller names have alternative STAR chimeric alignments. No specific PRAMEF rearrangement or protein target established.'),
10:('RNA_supported_nonproductive_orientation_unresolved','Exact raw junction support across unique25nt anchors, but3prime-3prime orientation and no caller peptide. Does not establish a functional CCND1 fusion, amplification, target or therapy.'),
11:('low_support_RNA_junction_unresolved','Five20+20/Q20 pairs and four25+25 pairs support this out-of-frame RNA junction; four strict50nt markers lie just1nt from a read end. A coding-to-lncRNA fusion/trans-splice/library chimera and tumor specificity remain unresolved.'),
12:('single_pair_RNA_junction_unresolved','Only one exact pair supports this alternate out-of-frame FOCAD-lncRNA junction. No independent molecular validation, tumor-specificity or confirmed neoantigen.'),
13:('single_pair_RNA_junction_unresolved','Only one exact pair supports this alternate out-of-frame FOCAD-lncRNA junction. No independent molecular validation, tumor-specificity or confirmed neoantigen.'),
14:('RNA_supported_nonproductive_orientation_unresolved','Raw junction corroborated, but FTH1-side paralog ambiguity and3prime-3prime orientation remain; no conventional fusion protein or specific DNA rearrangement established.'),
15:('repeat_or_paralog_confounded','Same marker as the other ACE partner calls; highly repetitive partner sequence prevents specific RPA3 assignment.3prime-3prime geometry is not a conventional coding fusion.'),
16:('repeat_or_paralog_confounded','Same marker as the other ACE partner calls; highly repetitive partner sequence prevents specific UMAD1 assignment.'),
17:('low_support_RNA_junction_unresolved','Five strict50nt/Q20 pairs/four full-pair sequence families support an out-of-frame LZTS2-CCND1 RNA junction. This is not evidence for a functional CCND1 activating fusion or a treatment indication.'),
18:('repeat_or_paralog_confounded','One anchor has tens of thousands of genome matches and all caller names have alternative STAR chimeric alignments.5prime-5prime orientation and stop-codon annotation do not establish a coding driver.'),
19:('repeat_or_paralog_confounded','Noncoding/repetitive upstream assignment and many genomic matches per anchor. Exact raw sequence does not resolve an intergenic driver or coding product.'),
20:('repeat_rich_RNA_junction_protein_effect_unresolved','Exact sequence is supported, with unique25nt genomic anchors but a highly GC/repeat-rich junction and short repeat-rich in-frame peptide. A genomic fusion versus repeat-length/placement/library explanation and functional protein are not established.'),
21:('paralog_or_overlapping_annotation_confounded','Same genomic breakpoints, read identifiers and exact marker as TUBB8P5-TUBB8P11; overlapping/pseudogene annotations and nonunique partner sequence. Not a distinct driver fusion.'),
22:('paralog_or_overlapping_annotation_confounded','Same genomic breakpoints, read identifiers and exact marker as TUBB8P5-FAM41C; overlapping/pseudogene annotations and nonunique partner sequence. No coding driver established.'),
23:('artifact_compatible_normal_reference_marker','The40nt junction matches ordinary ATXN3 transcripts and genomic sequence; repeat-side donor assignment is not specifically supported by this marker. Do not interpret its large raw count as fusion abundance.'),
24:('unresolved_ambiguous_repeat_consensus','Caller consensus contains unknown bases within a tandem-repeat context; a reference-only hypothesis yields one40nt/Q20 pair and no50nt pair. This is not a negative test for all alternate alleles/consensuses, and no tumor-specific ITD is established.'),
25:('repeat_or_paralog_confounded','Large counts arise from a poly-CAG marker with a nonunique TBP-side anchor and near-normal ATXN3 context; they are not specific evidence of a TBP-ATXN3 protein fusion. No validated coding product.'),
26:('artifact_compatible_normal_reference_marker','The40nt junction has a contiguous normal-genome match and a nonunique poly-CAG donor; same marker and breakpoint as the GLI4-labelled row.5prime-5prime/out-of-frame hypothesis, not an established EP400 fusion.'),
27:('artifact_compatible_normal_reference_marker','The40nt junction has a contiguous normal-genome match and a nonunique poly-CAG donor; same marker and breakpoint as the AC138696.1-labelled row.5prime-5prime/out-of-frame hypothesis, not an established EP400 fusion.')}
out=[]
for i,original in enumerate(p):
 r=dict(original);cid=r['candidate_id'];m=meta[cid]
 if cid in s:
  supplement=s[cid]
  # Replace only raw hypothesis evidence; retain actual Arriba source and support/name fields.
  for k,v in supplement.items():
   if k.startswith(('arm','junction_','left_anchor_','right_anchor_')) or k=='exact_pattern_arms_available':r[k]=v
  r['evidence_sequence_source']=sm[cid]['junction_hypothesis'];r['supplemental_source_names']=supplement['caller_named_pairs'];r['supplemental_source_names_recovered']=supplement['caller_named_pairs_recovered'];r['supplemental_name_definition']=sm[cid].get('name_source','Actual Arriba read_identifiers')
 else:r['evidence_sequence_source']='Caller-assembled consensus; exact simple one-breakpoint sequence';r['supplemental_source_names']='';r['supplemental_source_names_recovered']='';r['supplemental_name_definition']=''
 if i<28:status,reason=manual[i]
 else:
  if m['gene1']==m['gene2'] and 'duplication' in m['type']:
   status='discarded_self_backsplice_hypothesis_unresolved';reason='Reference-junction exact RNA support, where present, does not establish a genomic tandem duplication or activated kinase. Noncanonical/back-splice self-event was caller-filtered; no full-length isoform, protein product or functional domain retention validated.'
  elif m['gene1']==m['gene2'] and 'inversion' in m['type']:
   status='discarded_inverted_self_or_hairpin_hypothesis';reason='Exact RNA support, where present, has inverted same-gene breakends (3prime-3prime or5prime-5prime). Hairpin/library and nonproductive RNA explanations remain; this is not an established activating fusion or treatment target.'
  elif 'BAP1' in (m['gene1'],m['gene2']):
   status='discarded_BAP1_partner_not_second_hit';reason='Low-support, caller-filtered partner hypothesis; strict50nt support is0-4 pairs and repeat/sequence-family or short-anchor limits apply. No validated BAP1 structural second hit, allele phasing, protein product or biallelic inactivation is established.'
  else:
   status='discarded_partner_hypothesis_unresolved';reason='Caller-filtered candidate assessed with a joined-genome reference hypothesis, not an available assembled transcript. Zero exact support cannot exclude sequence variation, splice changes or non-template insertion; any positive support does not establish an oncogenic product.'
 r['source_name_definition']=r['supplemental_name_definition'] or 'Actual Arriba read_identifiers';r={k.replace('_caller_named_pairs','_source_named_pairs').replace('_caller_named_sequence_families','_source_named_sequence_families'):v for k,v in r.items()};r['review_status']=status;r['status_reason']=reason;r['research_only_not_clinical_validation']=True;r['matched_normal_available']=False;r['validated_actionability']='none established in this analysis';r['shared_marker_group']=';'.join(str(j+1) for j,g in enumerate(shared) if cid in g['candidates']);out.append(r)
# A marker-family count is sequence diversity, never a number of independent biological molecules.
fields=list(dict.fromkeys(k for r in out for k in r))
with (base/'candidate-reviewed.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(out)
(base/'candidate-reviewed.json').write_text(json.dumps(out,indent=2)+'\n')
counts=collections.Counter(r['review_status'] for r in out);(base/'review-status-counts.json').write_text(json.dumps(dict(counts),indent=2)+'\n');print(json.dumps(dict(counts),indent=2))
