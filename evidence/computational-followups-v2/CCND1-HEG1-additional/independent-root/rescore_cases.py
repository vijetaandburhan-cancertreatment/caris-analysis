from pathlib import Path
from dataclasses import asdict
import json,hashlib,time
from affine_window_dp import align
P=Path(__file__).resolve().parent.parent
all_cases=json.loads((P/'root-independent-rescore-cases.json').read_text())
# Before rescoring, deterministically take each base query/control's best exported
# proposed and best nonproposed case across both orientations. Validate winner
# scores, not an independently repeated exhaustive nomination/ranking.
groups={}
for c in all_cases:
 key=(c['query_id'][:-1],c['reference']['reference_kind']=='proposed')
 if key not in groups or (c['score'],c['case_id'])>(groups[key]['score'],groups[key]['case_id']):groups[key]=c
rows=[];started=time.time()
for key,c in sorted(groups.items()):
 result=align(c['query_sequence'],c['target_sequence']);assert result.score==c['score'],(c['case_id'],result.score,c['score'])
 rows.append({'query_base':key[0],'proposed':key[1],'case_id':c['case_id'],'reference':c['reference'],'query_length':len(c['query_sequence']),'window_length':len(c['target_sequence']),'primary_score':c['score'],'independent':asdict(result),'case_sha256':hashlib.sha256(json.dumps(c,sort_keys=True).encode()).hexdigest()})
 print(c['case_id'],result.score,flush=True)
out={'all_scores_match':True,'case_count':len(rows),'elapsed_seconds':round(time.time()-started,2),'selection_scope':'Independent rescore of best exported proposed and alternate cases per base query/control across orientations; no second genome nomination scan or exhaustive second ranking','input_sha256':hashlib.sha256((P/'root-independent-rescore-cases.json').read_bytes()).hexdigest(),'implementation_sha256':hashlib.sha256(Path(__file__).with_name('affine_window_dp.py').read_bytes()).hexdigest(),'results':rows}
Path(__file__).with_name('independent-patient-control-rescore.json').write_text(json.dumps(out,indent=2)+'\n')
print('DONE',len(rows),'scores matched in',out['elapsed_seconds'],'seconds')
