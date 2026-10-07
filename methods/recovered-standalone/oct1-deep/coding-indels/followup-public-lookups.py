"""Retrieve public catalog records after two primary API failures; no patient data uploaded."""
from pathlib import Path
import json,urllib.request,datetime
OUT=Path(__file__).resolve().parent/'public-catalog'
queries=[('APOL4-rs5845253-ncbi.json','https://api.ncbi.nlm.nih.gov/variation/v0/beta/refsnp/5845253'),('SLC27A3-gene-region-ensembl.json','https://rest.ensembl.org/overlap/region/human/1:153775407-153780157?feature=variation&content-type=application/json')]
for name,url in queries:
    try:
        with urllib.request.urlopen(urllib.request.Request(url,headers={'Accept':'application/json'}),timeout=50) as r:data=json.load(r)
        (OUT/name).write_text(json.dumps({'source':url,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'result':data},indent=2)+'\n')
        print(name,'success',len(data),flush=True)
    except Exception as e:print(name,'ERROR',repr(e),flush=True)
