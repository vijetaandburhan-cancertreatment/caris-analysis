"""Reproduce exact-allele official LOMO controls before any patient predictions."""
from pathlib import Path
import sys,json,csv,time,resource,hashlib,os
os.environ['CUDA_VISIBLE_DEVICES']=''
OUT=Path(__file__).resolve().parent;SRC=OUT/'public/code_and_dataset'
sys.path.insert(0,str(SRC))
import torch
torch.set_num_threads(2);torch.set_num_interop_threads(1)
from config_parser import Config
from model import Model
from seq_encoding import one_hot_PLUS_blosum_encode

def chains():
    rows={}
    for letter in ['A','B']:
        p=SRC/'dataset'/f'CLUATAL_OMEGA_{letter}_chains_aligned_FLATTEN.txt'
        rows.update(dict(r.split('\t')[:2] for r in p.read_text().splitlines()[1:]))
    return rows
CHAINS=chains()

def load_model(path):
    config=Config(str(SRC/'config_main.json'));config.device=torch.device('cpu')
    model=Model(config)
    state=torch.load(path,map_location='cpu',weights_only=True)
    model.load_state_dict(state,strict=True);model.eval()
    return model,config

def score(model,config,allele,peptide):
    assert 13<=len(peptide)<=25
    assert all(a in 'ACDEFGHIKLMNPQRSTVWY' for a in peptide)
    assert allele in ['DRB1*03:01','DRB1*15:02']
    arguments=[]
    for sequence,maxlen in [(CHAINS['DRA*01:01'],config.max_len_hla_A),(CHAINS[allele],config.max_len_hla_B),(peptide,config.max_len_pep)]:
        tensor,mask,length=one_hot_PLUS_blosum_encode(sequence,maxlen)
        arguments += [tensor.unsqueeze(0),mask.unsqueeze(0),torch.tensor([length])]
    with torch.inference_mode():value,attention=model(*arguments)
    value=float(value.item())
    assert 0<=value<=1
    assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss<2*1024**3
    return value

if __name__=='__main__':
    started=time.time();results=[];checks=[]
    for allele,folder in [('DRB1*03:01','DRA01:01-DRB103:01'),('DRB1*15:02','DRA01:01-DRB115:02')]:
        base=OUT/'public/Models/BD2016_LOMO'/folder
        model,config=load_model(base/'best_model.pytorch')
        data=list(csv.DictReader((base/'test_result.txt').open(),delimiter='\t'))
        unique={r['pep']:r for r in data if 13<=len(r['pep'])<=25}
        ordered=sorted(unique.values(),key=lambda r:(float(r['pred']),r['pep']))
        # Predeclared, deterministic spread over published score range.
        selected=ordered if len(ordered)<=25 else [ordered[round(i*(len(ordered)-1)/24)] for i in range(25)]
        vals=[]
        for r in selected:
            actual=score(model,config,allele,r['pep']);expected=float(r['pred']);delta=abs(actual-expected)
            vals.append(delta);checks.append({'allele':allele,'peptide':r['pep'],'expected_normalized_score':expected,'actual_normalized_score':actual,'absolute_difference':delta})
        passed=max(vals)<=1e-4
        results.append({'allele':allele,'model_scope':'Official BD2016 model trained with this target allele held out; this is not the BD2013 general model','unique_public_controls':len(selected),'max_absolute_score_difference':max(vals),'pass_tolerance':1e-4,'passed':passed})
        print(json.dumps(results[-1]),flush=True)
        del model
    result={'status':'PASS' if all(r['passed'] for r in results) else 'FAIL','controls':results,'checks':checks,'elapsed_seconds':time.time()-started,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'torch_version':torch.__version__,'score_direction':'Higher normalized score is stronger binding; source result_writer.py converts predicted_IC50_nM=50000**(1-score). Official main.py prints normalized score with an IC50 label; that print label must not be reused.','DRA_assumption':'Exact reference DRA*01:01 alpha sequence assumed, not measured here. Exact provisional DRB1*03:01 and DRB1*15:02 beta sequences used; no nearest-allele replacement.','limitations':'Public numerical regression is tool-execution validation, not a new accuracy benchmark or biological validation. No patient input was processed.'}
    (OUT/'public-controls.json').write_text(json.dumps(result,indent=2)+'\n')
    assert result['status']=='PASS','Stop: public control mismatch'
