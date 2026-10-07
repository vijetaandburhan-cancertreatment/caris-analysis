from pathlib import Path
import json,time,subprocess,sys,datetime
R=Path(__file__).resolve().parent;t=time.time()
while not (R/'patient-full-pass/run.json').exists():
 assert time.time()-t<1200,'Timed out waiting for full pass'
 time.sleep(5)
x=json.loads((R/'patient-full-pass/run.json').read_text());assert x['status']=='COMPLETE',x['status']
for name in ['count_compare.py','capture_original_names.py','compare_names.py']:
 print('START',name,datetime.datetime.now(datetime.timezone.utc).isoformat(),flush=True)
 with (R/(name+'.log')).open('w') as out:subprocess.run([sys.executable,str(R/name)],check=True,stdout=out,stderr=subprocess.STDOUT)
 print('COMPLETE',name,datetime.datetime.now(datetime.timezone.utc).isoformat(),flush=True)
