"""Public-only pinned DeepSeqPanII source/models; no patient input."""
from pathlib import Path
import urllib.request,urllib.parse,hashlib,json,time
OUT=Path(__file__).resolve().parent
repo='pcpLiu/DeepSeqPanII';commit='2161f5b8c873493c9d43c0201e88a7843d2999ec'
paths=['README.md','Models/README.md','code_and_dataset/config_main.json','code_and_dataset/config_parser.py','code_and_dataset/model.py','code_and_dataset/seq_encoding.py','code_and_dataset/main.py','code_and_dataset/result_writer.py','code_and_dataset/dataset/CLUATAL_OMEGA_A_chains_aligned_FLATTEN.txt','code_and_dataset/dataset/CLUATAL_OMEGA_B_chains_aligned_FLATTEN.txt','code_and_dataset/dataset/BD_2013_DATAPROVIDER_READY.txt','Models/benchmark_weekly/model_bd2013.pytorch','Models/benchmark_weekly/weekly_benchmark_results.txt']
for a in ['DRA01:01-DRB103:01','DRA01:01-DRB115:02']:
    paths += [f'Models/BD2016_LOMO/{a}/best_model.pytorch',f'Models/BD2016_LOMO/{a}/test_result.txt']
records=[]
for name in paths:
    url='https://raw.githubusercontent.com/'+repo+'/'+commit+'/'+urllib.parse.quote(name,safe='/')
    p=OUT/'public'/name;p.parent.mkdir(parents=True,exist_ok=True)
    if not p.exists():
        data=urllib.request.urlopen(url,timeout=60).read()
        assert len(data)<25_000_000
        p.write_bytes(data)
    records.append({'path':str(p.relative_to(OUT)),'url':url,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
assert sum(r['bytes'] for r in records)<55_000_000
(OUT/'public-source-manifest.json').write_text(json.dumps({'repo':repo,'commit':commit,'license_statement':'Official README License section says MIT; no account or eligibility click-through used.','downloaded_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'bytes':sum(r['bytes'] for r in records),'files':records},indent=2)+'\n')
print('Public-only files',len(records),'bytes',sum(r['bytes'] for r in records),flush=True)
