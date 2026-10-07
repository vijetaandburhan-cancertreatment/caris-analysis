from pathlib import Path
import json,collections
P=Path(__file__).resolve().parent.parent
queries=json.loads((P/'patient-queries.json').read_text());cases=json.loads((P/'independent-root/independent-patient-control-rescore.json').read_text())['results'];props={x['query_base']:x for x in cases if x['proposed']}
lookup={};spans={}
for q in queries:
 r=props[q['query_id']];a=r['independent'];coord=0;covered=[]
 for qb,rb in zip(a['query_aligned'],a['reference_aligned']):
  if rb!='-':
   if qb!='-':covered.append(coord)
   coord+=1
 assert coord==r['window_length']
 spans[q['query_id']]={'start0':min(covered),'end0':max(covered)+1,'orientation':r['case_id'].split('_')[0][-1],'score':a['score'],'case':r['case_id']}
 for m in q['members']:lookup[(m['pair_id'],m['mate'])]=q['query_id']
rows=[]
for p in json.loads((P/'pairs.json').read_text()):
 m=[spans[lookup[(p['pair_id'],k)]] for k in [1,2]]
 rows.append({'pair_id':p['pair_id'],'Q20':p['Q20_50nt_support'],'span':[min(x['start0'] for x in m),max(x['end0'] for x in m)],'opposite_orientations':m[0]['orientation']!=m[1]['orientation'],'mates':m})
q20=[x for x in rows if x['Q20']];patterns=collections.Counter(tuple(x['span']) for x in q20)
assert all(x['opposite_orientations'] for x in rows)
assert len(patterns)==6,patterns
out={'method':'Start/end intervals reconstructed from independent DP traceback, not owner coordinates. Full-reference gap-only ends excluded. Tie-selected placements can depend on alignment convention. Patterns are not proven original molecules.','Q20_start_end_patterns':{'%s-%s'%k:v for k,v in patterns.items()},'pattern_count':len(patterns),'all14_opposite_orientations':True,'pairs':rows}
(P/'independent-root/independent-pair-spans.json').write_text(json.dumps(out,indent=2)+'\n');print(out['Q20_start_end_patterns'])
