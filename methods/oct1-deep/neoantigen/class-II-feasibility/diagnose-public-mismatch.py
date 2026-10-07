"""Public-only one-example comparison with unchanged official main.run code."""
from pathlib import Path
import os,sys,json,importlib.util,io,contextlib,re
OUT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('controls',OUT/'run-public-controls.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
import torch
# File loading only: preserve required CPU and weights-only protections.
original_load=torch.load
def guarded_load(*args,**kwargs):
    kwargs['map_location']='cpu';kwargs['weights_only']=True
    return original_load(*args,**kwargs)
torch.load=guarded_load
def no_network(event,args):
    if event in ['socket.connect','socket.getaddrinfo']:raise RuntimeError('Network disabled in inference audit')
sys.addaudithook(no_network)
controls=json.loads((OUT/'public-controls.json').read_text())
record=max((r for r in controls['checks'] if r['allele']=='DRB1*15:02'),key=lambda r:r['absolute_difference'])
model_path=OUT/'public/Models/BD2016_LOMO/DRA01:01-DRB115:02/best_model.pytorch'
spec=importlib.util.spec_from_file_location('upstream_main',c.SRC/'main.py')
up=importlib.util.module_from_spec(spec);spec.loader.exec_module(up)
os.chdir(c.SRC)
buffer=io.StringIO()
with contextlib.redirect_stdout(buffer):
    up.run(str(model_path),'DRA*01:01',record['allele'],record['peptide'])
printed=buffer.getvalue()
value=float(printed.split(':')[-1].strip())
result={'scope':'One public example only; no patient input','allele':record['allele'],'peptide':record['peptide'],'released_expected_normalized_score':record['expected_normalized_score'],'wrapper_normalized_score':record['actual_normalized_score'],'unchanged_upstream_main_normalized_score':value,'upstream_and_wrapper_identical':abs(value-record['actual_normalized_score'])<1e-7,'official_main_print':printed.strip(),'only_interception':'torch.load forced map_location=cpu, weights_only=True; inference/model/encoding implementation unchanged','conclusion':'If upstream and wrapper agree, discrepancy predates wrapper; environment or released checkpoint/reference/control consistency remains unresolved. Do not run patient predictions.'}
(OUT/'public-mismatch-diagnosis.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
