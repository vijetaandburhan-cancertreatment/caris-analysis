import pathlib,json
P=pathlib.Path(__file__).parent/'whole-reference-alternatives'
r=json.loads((P/'scored-results.json').read_text());n=json.loads((P/'nominations.json').read_text());expected={'C_HEG1_parent':322,'C_CCND1_parent':322,'C_HEG1_three_substitutions':304,'C_HEG1_one_insertion':316,'C_HEG1_two_deletion':311,'C_MET_canonical_splice':322,'C_annotated_INS_IGF2_transcript':322}
rows=[]
for name,score in expected.items():
 x=next(x for x in r['results']if x['query']['oriented_id']==name+'+');truth=x['query']['truth'];nn=next(x for x in n['queries']if x['query']['oriented_id']==name+'+')
 hits=[w for w in nn['references'][truth['kind']]if w['reference_id']==truth['id'] and w['start0']<=truth['start0'] and w['end0']>=truth['end0']]
 best=max(w['score']for w in x['top_alternatives']);assert hits and best==score
 rows.append({'control':name,'expected_score':score,'best_alternative_score':best,'true_parent_window_nominated':True,'pass':True})
x=next(x for x in r['results']if x['query']['oriented_id']=='C_synthetic_proposed_join+');assert max(m['score']for m in x['proposed_models'])==322 and max(m['score']for m in x['top_alternatives'])==148
rows.append({'control':'C_synthetic_proposed_join','proposed_score':322,'best_normal_score':148,'pass':True})
x=next(x for x in r['results']if x['query']['oriented_id']=='C_repeat_nonunique+');nn=next(x for x in n['queries']if x['query']['oriented_id']=='C_repeat_nonunique+')
assert len(nn['references']['genome'])>1 and len(nn['references']['transcript'])>1
rows.append({'control':'C_repeat_nonunique','genome_windows':len(nn['references']['genome']),'transcript_windows':len(nn['references']['transcript']),'best_normal_score':max(m['score']for m in x['top_alternatives']),'pass':'many candidate placements; no claim of a unique origin or exact 161nt repeat match'})
(P/'control-validation.json').write_text(json.dumps({'controls':rows,'seed_caps':n['seed_caps'],'window_caps':n['window_caps'],'all_planned_tests_pass':True},indent=2)+'\n')
print(json.dumps(rows))
