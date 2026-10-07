#!/usr/bin/env python3
import collections,gzip,hashlib,importlib.util,json,pathlib,sys,time
base=pathlib.Path(__file__).resolve().parent.parent
spec=importlib.util.spec_from_file_location('declared_triage',base/'read-level-audit/inventory_candidates_missing_coverage.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
run=base/'patient-D8';out=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct4-fusion-read-audit/patient-inventory');out.mkdir(parents=True,exist_ok=True)
selected=[];counts=collections.Counter();confidence=collections.Counter();t=time.time();sources={}
with gzip.open(out/'all-triaged-candidates.jsonl.gz','wt',compresslevel=3) as f:
 for filename,kind in [('fusions.tsv','accepted'),('fusions.discarded.tsv','discarded')]:
  source=run/filename;h=hashlib.sha256()
  with open(source,'rb') as b:
   for chunk in iter(lambda:b.read(1024*1024),b''):h.update(chunk)
  sources[filename]={'bytes':source.stat().st_size,'sha256':h.hexdigest()}
  for row in mod.read(source,kind):
   counts[kind]+=1;f.write(json.dumps(row,separators=(',',':'))+'\n')
   if kind=='accepted':confidence[row['confidence']]+=1
   if row['review_reasons']:selected.append(row)
(out/'review-candidates.json').write_text(json.dumps(selected,indent=2)+'\n')
summary={'counts':dict(counts),'accepted_confidence':dict(confidence),'selected':len(selected),'selected_by_source':dict(collections.Counter(r['source_set'] for r in selected)),'selected_with_junction_patterns':sum(bool(r['junction_patterns']) for r in selected),'unique_caller_names':len(set(n for r in selected for n in r['read_names'])),'sources':sources,'selection_script':str(base/'inventory_candidates.py'),'selection_script_sha256':hashlib.sha256((base/'inventory_candidates.py').read_bytes()).hexdigest(),'adaptation':'Declared read() generator with a single ingestion correction: missing coverage1/coverage2 dot retained as null. Stream complete inventory as gzip JSONL. Identical selection predicates.','elapsed_seconds':time.time()-t}
(out/'inventory.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
