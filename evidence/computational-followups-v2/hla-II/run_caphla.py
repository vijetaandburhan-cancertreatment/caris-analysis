#!/usr/bin/env python3
"""Pinned upstream CapHLA models/encoding; CPU wrapper removes multiprocessing only.
No input/output network operations; refuse socket activity during inference.
Public regression predeclared normalized-score max absolute tolerance: 1e-5.
"""
import os,sys,pathlib,json,hashlib,time,resource,argparse
os.environ['OMP_NUM_THREADS']='2';os.environ['MKL_NUM_THREADS']='2'
import numpy as np,pandas as pd,torch
from torch.utils.data import DataLoader
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'CapHLA'))
from BA_model import CapHLA_BA
from EL_model import CapHLA_EL
from utils import _Pep_MHC_dataset,aa_dict_one_hot,predict_ba,predict_ms

def deny_network(event,args):
 if event.startswith('socket.'): raise RuntimeError('Network prohibited during model inference: '+event)

def predict(df,modeldir=ROOT/'CapHLA',batch_size=32,library=None):
 ldf=pd.read_csv(library or modeldir/'HLA_library.csv'); lib=ldf.set_index('Allele Name')['MHC pseudo-seq'].to_dict()
 for s in df.peptide:
  assert 7<=len(s)<=25 and set(s)<=set('ACDEFGHIKLMNPQRSTVWY')
 h=[lib[a] for a in df['Allele Name']]; assert all(len(x)==34 for x in h)
 dataset=_Pep_MHC_dataset([[list(s+'X'*(25-len(s))) for s in df.peptide],[list(x) for x in h]],aa_dict_one_hot)
 iterator=DataLoader(dataset,batch_size=batch_size,shuffle=False,num_workers=0)
 ans=df.copy(); foldout={};torch.set_num_threads(2)
 for name,cls,fun in [('el',CapHLA_EL,predict_ms),('ba',CapHLA_BA,predict_ba)]:
  folds=[]
  for fold in range(5):
   net=cls().to('cpu');net.load_state_dict(torch.load(modeldir/'params'/f'{name}_fold{fold}.params',map_location='cpu',weights_only=True),strict=True);net.eval()
   values=fun(net,iterator,torch.device('cpu'));folds.append(values);ans[f'{name}_fold{fold}']=values
  ans['presentation_score' if name=='el' else 'affinity_score']=np.asarray(folds).mean(axis=0)
 return ans

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--output',required=True);ap.add_argument('--expected');ap.add_argument('--model-dir');ap.add_argument('--library');ap.add_argument('--batch-size',type=int,default=32);a=ap.parse_args()
 sys.addaudithook(deny_network);start=time.time();df=pd.read_csv(a.input);out=predict(df,pathlib.Path(a.model_dir) if a.model_dir else ROOT/'CapHLA',a.batch_size,a.library);out.to_csv(a.output,index=False)
 meta={'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'python':sys.version,'torch':torch.__version__,'numpy':np.__version__,'pandas':pd.__version__,'rows':len(out),'elapsed_seconds':time.time()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'batch_size':a.batch_size,'library':str(a.library),'network_blocked':True,'model_dir':str(a.model_dir or ROOT/'CapHLA')}
 if a.expected:
  expected=pd.read_csv(a.expected);assert len(expected)==len(out);assert list(zip(expected.peptide,expected['Allele Name']))==list(zip(out.peptide,out['Allele Name']))
  meta['predeclared_max_absolute_tolerance']=1e-5; meta['controls']={}
  for c in ['presentation_score','affinity_score']:
   dif=abs(out[c]-expected[c]);meta['controls'][c]={'max_absolute_error':float(dif.max()),'rmse':float(np.sqrt((dif**2).mean())),'passed':bool(dif.max()<=1e-5)}
  meta['passed']=all(x['passed'] for x in meta['controls'].values())
 pathlib.Path(a.output+'.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta,indent=2))
