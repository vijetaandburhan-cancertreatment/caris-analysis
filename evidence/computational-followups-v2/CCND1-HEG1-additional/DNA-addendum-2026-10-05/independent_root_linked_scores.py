from pathlib import Path
import sys,json,hashlib,datetime
from dataclasses import asdict
W=Path(__file__).resolve().parent
sys.path.insert(0,str(W.parent/'independent-root'))
from affine_window_dp import align
src=Path(sys.argv[1]) if len(sys.argv)>1 else W/'linked-root-rescore-cases.json'
rows=[]
for c in json.loads(src.read_text()):
 r=align(c['query_sequence'],c['target_sequence'])
 assert r.score==c['score'],(c['case_id'],r.score,c['score'])
 rows.append({'case_id':c['case_id'],'score':r.score,'independent_alignment':asdict(r),'reference':c['reference']})
out={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'case_count':len(rows),'all_scores_match':True,'scope':'Independent standard-library dynamic-programming rescore of supplied windows; not an independent genome nomination search, and scores are not probabilities.','rows':rows}
dst=W/('independent-'+src.stem+'.json');dst.write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'output':str(dst),'case_count':len(rows),'all_scores_match':True,'scores':{x['case_id']:x['score'] for x in rows}}))
