import pathlib,json,pandas as pd,numpy as np,hashlib,datetime
P=pathlib.Path(__file__).resolve().parent
manifest=json.loads((P/'patient-peptide-manifest.json').read_text());wins={w['window_id']:w for w in manifest['windows']};ref=json.loads((P/'normal-reference-exact-check.json').read_text());frames={v:pd.read_csv(P/f'patient-{v}-scores.csv') for v in ['v2','v1']}
key=['window_id','Allele Name'];pairs=None
for v,d in frames.items():
 assert len(d)==manifest['input_rows'];assert d[key+['sequence_type']].duplicated().sum()==0
 m=d[d.sequence_type=='mutant'].set_index(key);w=d[d.sequence_type=='wt'].set_index(key);assert set(m.index)==set(w.index)
 r=m[['peptide','gene','length','aa_start1','aa_end1','altered_positions_in_peptide1','reaches_predicted_stop','HLA_caveat']].copy();r['WT_peptide']=w.peptide
 for kind,sub in [('mutant',m),('wt',w)]:
  for score,pre in [('presentation_score','el'),('affinity_score','ba')]:
   r[f'{v}_{kind}_{pre}']=sub[score]
   vals=sub[[f'{pre}_fold{i}' for i in range(5)]]
   r[f'{v}_{kind}_{pre}_fold_min']=vals.min(axis=1);r[f'{v}_{kind}_{pre}_fold_max']=vals.max(axis=1)
 for typ in ['el','ba']:r[f'{v}_mutant_minus_wt_{typ}']=r[f'{v}_mutant_{typ}']-r[f'{v}_wt_{typ}']
 if pairs is None:pairs=r
 else:pairs=pairs.join(r[[c for c in r if c.startswith(v+'_')]])
pairs=pairs.reset_index();rows=[]
for r in pairs.to_dict('records'):
 w=wins[r['window_id']];nov=[];shared=[]
 for i in range(len(w['mutant'])-8):
  s=w['mutant'][i:i+9]
  (nov if any(i+1<=q<=i+9 for q in w['altered_positions_in_peptide1']) else shared).append(s)
 r.update({'mutation_containing_9mer_windows':len(nov),'wholly_WT_9mer_windows':len(shared),'mutation_containing_9mer_sequences':';'.join(nov),'mutated_residues_in_peptide':len(w['altered_positions_in_peptide1']),'official_HLA_pseudosequence_has_X':next(a for a in manifest['alleles'] if a['name']==r['Allele Name'])['pseudo_sequence'].count('X')>0})
 for v in ['v1','v2']:r[f'{v}_both_scores_above_author_0_5_cutpoints']=bool(r[f'{v}_mutant_el']>.5 and r[f'{v}_mutant_ba']>.5)
 rows.append(r)
d=pd.DataFrame(rows);d.to_csv(P/'patient-paired-comparisons.csv',index=False)
# Always show strongest EL result per gene/HLA and strongest BA result, not only favorable WT contrast.
by=[]
for (gene,a),g in d.groupby(['gene','Allele Name']):
 for criterion,col in [('highest_v2_EL','v2_mutant_el'),('highest_v2_BA','v2_mutant_ba'),('largest_v2_EL_mutant_minus_WT','v2_mutant_minus_wt_el')]:
  r=g.sort_values([col,'window_id'],ascending=[False,True]).iloc[0].to_dict();r['selection_criterion']=criterion;by.append(r)
pd.DataFrame(by).to_csv(P/'per-HLA-candidate-comparison.csv',index=False)
summary=[]
for gene,g in d.groupby('gene'):
 row={'gene':gene,'distinct_mutant_peptides':int(g.window_id.nunique()),'HLA_combinations':int(g['Allele Name'].nunique()),'peptide_HLA_pairs':len(g),'max_v2_EL':float(g.v2_mutant_el.max()),'max_v2_BA':float(g.v2_mutant_ba.max()),'mutant_v2_EL_gt_WT':int((g.v2_mutant_minus_wt_el>0).sum()),'mutant_v2_BA_gt_WT':int((g.v2_mutant_minus_wt_ba>0).sum()),'both_author_cutpoints_v2':int(g.v2_both_scores_above_author_0_5_cutpoints.sum()),'both_author_cutpoints_v1':int(g.v1_both_scores_above_author_0_5_cutpoints.sum()),'both_author_cutpoints_both_releases':int((g.v2_both_scores_above_author_0_5_cutpoints&g.v1_both_scores_above_author_0_5_cutpoints).sum()),'max_EL_difference_between_releases':float(abs(g.v2_mutant_el-g.v1_mutant_el).max()),'max_BA_difference_between_releases':float(abs(g.v2_mutant_ba-g.v1_mutant_ba).max())}
 summary.append(row)
out={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'model':'CapHLA; v2 primary, v1 release sensitivity; same algorithm family, not independent method validation','gene_summary':summary,'peptide_input_manifest_sha256':hashlib.sha256((P/'patient-peptide-manifest.json').read_bytes()).hexdigest(),'input_rows':manifest['input_rows'],'comparison_rows':len(d),'top_rows':by,'score_semantics':'EL = model softmax output trained on eluted-ligand/negative examples, not calibrated probability of patient tumor presentation. BA = raw normalized affinity-regression score, higher means predicted stronger binding; neither nM nor percentile rank. Five model folds averaged as original code. No probability of vaccine or treatment success.','normal_reference':{'release':'GENCODE50','protein_records':ref['reference_records_scanned'],'long_exact_matches':len(ref['long_peptide_exact_matches']),'mutation_containing9mers':ref['novel_sequence_core_count'],'mutation_containing9mers_with_exact_reference_match':sum(x['reference_occurrences']>0 for x in ref['cores'].values() if x['mutation_containing'])}}
(P/'summary.json').write_text(json.dumps(out,indent=2));print(json.dumps(summary,indent=2));print(pd.DataFrame(by).query("selection_criterion=='highest_v2_EL'")[['gene','Allele Name','peptide','WT_peptide','v2_mutant_el','v2_wt_el','v2_mutant_ba','v2_wt_ba','v1_mutant_el','v1_mutant_ba','wholly_WT_9mer_windows']].to_string(index=False))
