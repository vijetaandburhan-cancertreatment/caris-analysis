"""Local affinity comparison for selecting experimental controls, not safety prediction."""
from pathlib import Path
import json,os,time,csv
os.environ['CUDA_VISIBLE_DEVICES']=''
os.environ['MHCFLURRY_DEVICE']='cpu'
import torch
torch.set_num_threads(2)
from mhcflurry import Class1AffinityPredictor
OUT=Path(__file__).resolve().parent
BASE=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/mhcflurry-models')
data=json.loads((OUT/'bap1-normal-sequence-neighbors.json').read_text())
peptides=['IEERKGLYL','EERKGLYL','IEERRSRRK','EERRSRRK']+[x['normal_reference_peptide'] for x in data['results']]
assert len(peptides)==len(set(peptides))
rows=[]
for release in ['2.3.0','2.2.0']:
    model=BASE/release/'models_class1_presentation/models/affinity_predictor'
    p=Class1AffinityPredictor.load(str(model))
    vals=p.predict(peptides,allele='HLA-B*50:01')
    ranks=p.percentile_ranks(vals,allele='HLA-B*50:01')
    for peptide,value,rank in zip(peptides,vals,ranks):
        rows.append({'peptide':peptide,'allele':'HLA-B*50:01','model_release':release,'predicted_affinity_nM':float(value),'affinity_percentile':float(rank),'type':'BAP1 altered' if peptide in peptides[:2] else 'BAP1 normal counterpart' if peptide in peptides[2:4] else 'normal-reference two-substitution neighbor'})
    del p
result={'scope':'Model-only affinity comparison with 24 normal-reference equal-length two-substitution neighbors. No cross-reactivity, safety, peptide presentation, immunogenicity or clinical efficacy conclusion. Neighbors are potential experimental control peptides, not predicted off-targets. Reference proteome excludes population variation/noncanonical translation and only equal-length Hamming <=2 neighbors were searched.','rows':rows}
(OUT/'bap1-normal-neighbor-affinity.json').write_text(json.dumps(result,indent=2)+'\n')
with (OUT/'bap1-normal-neighbor-affinity.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
for r in sorted((x for x in rows if x['model_release']=='2.3.0'),key=lambda x:x['predicted_affinity_nM']):print(r,flush=True)
