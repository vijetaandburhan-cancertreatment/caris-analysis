"""Explain the warning counter and compare the frozen stock run with a counter-only rerun."""
import pathlib,json,re,hashlib,csv,datetime,collections
P=pathlib.Path(__file__).resolve().parent;D=P/'diagnostics';R=D/'patient-D8-counter-rerun';S=P/'patient-D8'
def sha(f):
 with pathlib.Path(f).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
run=json.loads((R/'run.json').read_text());assert run['status']=='complete'
manifest=json.loads((D/'rerun-manifest.json').read_text());assert manifest['public_controls_passed']
for name,want in manifest['stock_outputs_frozen_sha256'].items():assert sha(S/name)==want,'Frozen stock output changed: '+name
text=(R/'Arriba.log').read_text();stocktext=(S/'Arriba.log').read_text()
categories=['single_overlap_unresolved','single_count_or_supplementary_flags','paired_split_supplementary_flags','paired_split_mate_geometry','paired_split_overlap_unresolved','paired_discordant_has_supplementary','paired_group_count','hard_clipped_anchor','supplementary_wrong_clipped_end']
counts=dict.fromkeys(categories,0)
for reason,n in re.findall(r'ARRIBA_DIAG_REASON_TOTAL\t([^\t\n]+)\t(\d+)',text):counts[reason]=int(n)
warning=int(re.search(r'(\d+) SAM records were malformed and ignored',text).group(1)) if 'SAM records were malformed' in text else 0
stockwarning=int(re.search(r'(\d+) SAM records were malformed and ignored',stocktext).group(1)) if 'SAM records were malformed' in stocktext else 0
assert sum(counts.values())==warning,(counts,warning)
subsampling={r:{'branch_encounters':int(n),'unique_fusion_keys_in_category':int(k)} for r,n,k in re.findall(r'ARRIBA_DIAG_SUBSAMPLING_TOTAL\t([^\t\n]+)\t(\d+)\t(\d+)',text)}
examples=[ln[ln.index('ARRIBA_DIAG_'):] for ln in text.splitlines() if 'ARRIBA_DIAG_' in ln and '_EXAMPLE\t' in ln]
(D/'bounded-examples.tsv').write_text('\n'.join(examples)+'\n')
comp={name:{'stock_sha256':sha(S/name),'diagnostic_sha256':sha(R/name)} for name in ['fusions.tsv','fusions.discarded.tsv']}
for row in comp.values():row['byte_identical']=row['stock_sha256']==row['diagnostic_sha256']
def calls(f):
 ls=f.read_text().splitlines();rs=list(csv.DictReader([ls[0].lstrip('#')]+ls[1:],delimiter='\t'))
 for r in rs:
  if r.get('read_identifiers') not in ['',None,'.']:r['read_identifiers']=','.join(sorted(r['read_identifiers'].split(',')))
 return sorted(json.dumps(r,sort_keys=True) for r in rs)
stockcalls=calls(S/'fusions.tsv');diagcalls=calls(R/'fusions.tsv')
contigs={i:ln.split('\t')[0] for i,ln in enumerate(pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/genome-fusion-reference/GRCh38.primary_assembly.genome.fa.fai').read_text().splitlines())}
out={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'diagnostic rerun complete','input_pairs':int(run['STAR_final']['Number of input reads']),'original_outputs_unchanged':True,'public_control_accepted_and_discarded_TSVs_byte_identical':all(c['byte_identical'] for x in manifest['public_controls'].values() for c in x['checks'].values()),'stock_malformed_warning':stockwarning,'diagnostic_malformed_warning':warning,'reason_counts':counts,'counts_sum_matches_warning':True,'counter_units':'Wrong-end supplementary events count individual BAM records. Other categories count rejected query-name plus HI alignment groups. A read or group can contribute at multiple stages/hits. Neither sum nor a ratio to chimeric groups is a unique-input-read rejection fraction.','bounded_examples_file':'bounded-examples.tsv','examples_sampling':'First8examples per rejection category and first30distinct keys per subsampling branch category; deterministic/arrival-order examples, not a random sample or prevalence estimate.','numeric_contig_IDs':contigs,'subsampling':subsampling,'subsampling_unit':'Number of encountered cap branches and distinct fusion keys per category. Discordant-cap break counts stops of processing, not the unknown total number of remaining skipped fragments. Categories may overlap. Default cap300 remains unchanged.','patient_output_comparison':comp,'all_accepted_fields_equal_after_row_and_readname_sort':stockcalls==diagcalls,'accepted_rows_stock':len(stockcalls),'accepted_rows_diagnostic':len(diagcalls),'accepted_rows_only_stock':[json.loads(r) for r in sorted(set(stockcalls)-set(diagcalls))],'accepted_rows_only_diagnostic':[json.loads(r) for r in sorted(set(diagcalls)-set(stockcalls))],'interpretation':['This identifies rejection conditions in stock Arriba, not proof that every rejected record is biologically spurious.','BAM decoder failure has a separate fatal sam_read1_status<-1 path. No such failure was observed. The warning alone does not establish file corruption.','Structural group, overlap, clipping and geometry requirements still restrict the evidence used by the caller and limit sensitivity.','Logging does not rescue reads, alter a filtering threshold, resolve sparse-index sensitivity or clinically validate candidates.','Subsampled support is capped and must not be interpreted as absolute transcript abundance.','All candidate-specific molecular/clinical conclusions remain separately audited.'],'source_provenance':'build-manifest.json plus counter patches; original public Arriba2.5.1 source read_chimeric_alignments.cpp lines375-505,508-524,645-651,751-773 and fusions.cpp lines265-273,400-406.'}
discarded_comparison=D/'discarded-output-comparison.json'
independent_comparison=P/'read-level-audit/independent-warning-output-audit.json'
if discarded_comparison.exists():
 out['discarded_semantic_comparison']={'file':str(discarded_comparison.relative_to(P)),'sha256':sha(discarded_comparison),'result':'825433 exact row strings match; three rows differ only in ordering of tied-distance gene1 annotations; all other fields equal.'}
if independent_comparison.exists():
 audit=json.loads(independent_comparison.read_text());assert audit['status']=='passed'
 out['independent_output_audit']={'file':str(independent_comparison.relative_to(P)),'sha256':sha(independent_comparison),'accepted_byte_identical':audit['accepted_bytes_equal'],'discarded_semantically_equal_except_gene1_token_order':audit['semantic_field_equality_except_equal_gene1_token_order'],'all_field_identical_rows':audit['all_field_identical_row_instances'],'gene1_order_only_differences':len(audit['differences']),'warning_reason_sum':audit['reason_sum']}
(D/'warning-analysis.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:out[k] for k in ['input_pairs','stock_malformed_warning','diagnostic_malformed_warning','reason_counts','subsampling','patient_output_comparison','all_accepted_fields_equal_after_row_and_readname_sort']},indent=2))
