"""Assert and export the compact benchmark conclusions from completed source-backed results."""
from pathlib import Path
import json,hashlib,collections,csv,datetime
P=Path(__file__).resolve().parent
r=json.loads((P/'results.json').read_text());d=json.loads((P/'design-manifest.json').read_text())
assert r['runs_completed']==r['expected_runs']==45
runs={x['run']:x for x in r['runs']}
comp=r['paired_compact_comparisons'];assert len(comp)==20
assert all(x['accepted_D1']==x['accepted_D8'] for x in comp)
changed=[x for x in comp if x['Arriba_names_only_D1'] or x['Arriba_names_only_D8']]
assert len(changed)==4 and all(x['Arriba_names_only_D1']==['BCR-ABL1-28'] and x['Arriba_names_only_D8']==[] for x in changed)
ids=collections.defaultdict(list)
for x in d['datasets']:ids[tuple(z['uncompressed_sha256'] for z in x['files'])].append(x['id'])
assert len(ids)==19
trim=[]
for setting in ['compactD1','compactD8','fullD8']:
    base=runs[setting+'__strong_original'];caps=[12] if setting=='fullD8' else [12,20,30]
    assert base['accepted_support_count_names']==3
    for cap in caps:
        t=runs[setting+f'__strong_trim_anchor{cap}'];assert t['accepted_support_count_names']==(3 if cap==30 else 2)
        lost=sorted(set(base['accepted_support_names'])-set(t['accepted_support_names']))
        assert lost==([] if cap==30 else ['BCR-ABL1-12'])
        trim.append({'setting':setting,'condition_individual_read_anchor_cap':cap,'untrimmed_support_names':base['accepted_support_names'],'trimmed_support_names':t['accepted_support_names'],'lost_support_names':lost,'interpretation':'Paired-read trimming condition. Complementary mates can merge and restore long effective anchors; no universal anchor threshold inferred.'})
a=runs['compactD1__support_all'];b=runs['compactD8__support_all'];c=runs['fullD8__support_all']
assert [x['input_pairs'] for x in [a,b,c]]==[45,45,45]
assert [x['accepted_support_count_names'] for x in [a,b,c]]==[10,9,10]
lost=sorted(set(a['accepted_support_names'])-set(c['accepted_support_names']));gained=sorted(set(c['accepted_support_names'])-set(a['accepted_support_names']))
assert lost==['BCR-ABL1-28'] and gained==['BCR-ABL1-18']
n=runs['compactD1__no_known_recovery__support_01_selection_1'];y=runs['compactD1__support_01_selection_1']
assert not n['target_accepted'] and not y['target_accepted'] and n['target_discarded_filters']==y['target_discarded_filters']==['in_vitro']
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
out={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'design_manifest_sha256':sha(P/'design-manifest.json'),'results_sha256':sha(P/'results.json'),'planned_cases_completed':45,'ordinary_caller_table_cases':sum(x['caller_outcome']=='completed_with_call_tables' for x in runs.values()),'explicit_no_chimeric_input_exit1_cases':sum(x['caller_outcome'].startswith('no_chimeric_input') for x in runs.values()),'compact_matched_condition_labels':20,'unique_paired_input_sets':19,'duplicate_input_labels':[v for v in ids.values() if len(v)>1],'compact_accepted_decision_agrees_in_all_conditions':True,'sparseD8_support_loss_cases':changed,'paired_read_trimming':trim,'all13_source_pairs_comparison':{'input_pairs_each':45,'normal_parent_pairs':32,'source_pairs':13,'compactD1_support_names':a['accepted_support_names'],'compactD8_support_names':b['accepted_support_names'],'fullD8_support_names':c['accepted_support_names'],'fullD8_vs_compactD1_lost':lost,'fullD8_vs_compactD1_gained':gained},'no_known_recovery_check':{'source_pair':'BCR-ABL1-4','with_k':y['target_discarded_filters'],'without_k':n['target_discarded_filters'],'scope':'One predeclared case; candidate-level filter, no general conclusion about importance of known-fusion priors.'},'main_limitations':['Twenty condition labels are nineteen unique inputs, not independent biological replicates.','Small32-parent background cannot calibrate confidence/e-values, patient-library abundance, clinical sensitivity, specificity or detection limits.','Compact reference omits competing loci; no full-reference D1/D8 comparison was performed.','Paired-read trimming can leave long effective merged-fragment anchors; the nominal cap is not a validated fragment-level threshold.','No patient rerun, clinical reclassification or new treatment target is established.']}
assert out['ordinary_caller_table_cases']==38 and out['explicit_no_chimeric_input_exit1_cases']==7
(P/'summary.json').write_text(json.dumps(out,indent=2))
with (P/'summary.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t');w.writerow(['condition','D1_accepted','D8_accepted','D1_support_names','D8_support_names','D1_only','D8_only'])
    for x in comp:
        name=x['dataset'];ra=runs['compactD1__'+name];rb=runs['compactD8__'+name]
        w.writerow([name,x['accepted_D1'],x['accepted_D8'],','.join(ra['accepted_support_names']),','.join(rb['accepted_support_names']),','.join(x['Arriba_names_only_D1']),','.join(x['Arriba_names_only_D8'])])
print(json.dumps({k:out[k] for k in ['planned_cases_completed','ordinary_caller_table_cases','explicit_no_chimeric_input_exit1_cases','compact_matched_condition_labels','unique_paired_input_sets','compact_accepted_decision_agrees_in_all_conditions']},indent=2))
