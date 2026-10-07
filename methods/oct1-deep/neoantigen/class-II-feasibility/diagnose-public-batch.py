"""Public-only check of upstream training batch size 64 versus main batch size1."""
from pathlib import Path
import sys,json,importlib.util,resource,time
OUT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('controls',OUT/'run-public-controls.py')
c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
import torch
def no_network(event,args):
    if event in ['socket.connect','socket.getaddrinfo']:raise RuntimeError('Network disabled')
sys.addaudithook(no_network)
r=json.loads((OUT/'public-mismatch-diagnosis.json').read_text())
config=c.Config(str(c.SRC/'config_main.json'));config.device=torch.device('cpu');config.config['Training']['batch_size']=64
model=c.Model(config)
model.load_state_dict(torch.load(OUT/'public/Models/BD2016_LOMO/DRA01:01-DRB115:02/best_model.pytorch',map_location='cpu',weights_only=True),strict=True)
model.eval();arguments=[]
for sequence,maxlen in [(c.CHAINS['DRA*01:01'],274),(c.CHAINS['DRB1*15:02'],291),(r['peptide'],25)]:
    tensor,mask,length=c.one_hot_PLUS_blosum_encode(sequence,maxlen)
    arguments += [tensor.unsqueeze(0).repeat(64,1,1),mask.unsqueeze(0).repeat(64,1),torch.tensor([length]*64)]
started=time.time()
with torch.inference_mode():score,attention=model(*arguments)
values=[float(x) for x in score.flatten()]
res={'scope':'Repeated one public control peptide only; no patient data','batch_size':64,'batch1_score':r['wrapper_normalized_score'],'batch64_score_min':min(values),'batch64_score_max':max(values),'expected_released_score':r['released_expected_normalized_score'],'maximum_batch_effect':max(abs(x-r['wrapper_normalized_score']) for x in values),'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'elapsed_seconds':time.time()-started}
(OUT/'public-batch-diagnosis.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
