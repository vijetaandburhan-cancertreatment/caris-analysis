"""Same local HLA-I models, updated public reference neighbors. Not a safety test."""
from pathlib import Path
import csv, json, os, sys, math, time, resource, hashlib
from datetime import datetime, timezone
os.environ['CUDA_VISIBLE_DEVICES'] = ''
os.environ['MHCFLURRY_DEVICE'] = 'cpu'
import torch
torch.set_num_threads(2)
from mhcflurry import Class1AffinityPredictor

def no_network(event, args):
    if event in ('socket.connect', 'socket.connect_ex', 'urllib.Request'):
        raise RuntimeError('Network disabled for local prediction')
sys.addaudithook(no_network)
OUT = Path(__file__).resolve().parent
MODELS = Path('/Users/burhanazeem/.local/share/codex/caris-analysis/mhcflurry-models')
start = time.monotonic()
neighbors = list(csv.DictReader((OUT/'gencode50-neighbors.tsv').open(),delimiter='\t'))
seq = [r['reference_peptide'] for r in neighbors]
assert len(seq)==31 and len(set(seq))==31
assert all(r['query']=='EERKGLYL' and int(r['substitutions'])==2 for r in neighbors)
controls = ['IEERKGLYL','EERKGLYL','IEERRSRRK','EERRSRRK']
peptides = controls + seq
old = json.loads((OUT.parent/'bap1-normal-neighbor-affinity.json').read_text())
lookup = {(r['model_release'],r['peptide']):r for r in old['rows']}
rows=[]; regression=[]
for release in ['2.3.0','2.2.0']:
    p=Class1AffinityPredictor.load(str(MODELS/release/'models_class1_presentation/models/affinity_predictor'))
    vals=p.predict(peptides,allele='HLA-B*50:01')
    ranks=p.percentile_ranks(vals,allele='HLA-B*50:01')
    for peptide,value,rank in zip(peptides,vals,ranks):
        row={'peptide':peptide,'allele':'HLA-B*50:01','model_release':release,
             'predicted_affinity_nM':float(value),'affinity_percentile':float(rank),
             'role':'BAP1 altered' if peptide in controls[:2] else 'BAP1 normal counterpart' if peptide in controls[2:] else 'reference two-substitution control',
             'in_prior_GENCODE37_run':(release,peptide) in lookup}
        rows.append(row)
        if (release,peptide) in lookup:
            previous=lookup[(release,peptide)]['predicted_affinity_nM']
            passed=math.isclose(float(value),previous,rel_tol=1e-4,abs_tol=0.1)
            regression.append({'model':release,'peptide':peptide,'old_nM':previous,'new_nM':float(value),'pass':passed})
    del p
passed=all(r['pass'] for r in regression)
data={'created_utc':datetime.now(timezone.utc).isoformat(),
      'scope':'Local model binding predictions for reference-sequence control selection only. No cross-reactivity, safety, presentation or efficacy conclusion.',
      'source_neighbors_sha256':hashlib.sha256((OUT/'gencode50-neighbors.tsv').read_bytes()).hexdigest(),
      'same_model_replay_tolerance':{'relative':1e-4,'absolute_nM':0.1},
      'prior_GENCODE37_replay_pass':passed,'prior_replay':regression,'rows':rows,
      'elapsed_seconds':round(time.monotonic()-start,2),'peak_RSS_bytes_macOS':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
(OUT/'reference-control-affinity.json').write_text(json.dumps(data,indent=2)+'\n')
with (OUT/'reference-control-affinity.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps({'passed':passed,'replayed':len(regression),'new_peptides':7,'rows':len(rows),'seconds':data['elapsed_seconds'],'RSS_bytes':data['peak_RSS_bytes_macOS']}))
if not passed:raise RuntimeError('Same-model replay failed; do not interpret updated predictions')
print(json.dumps([r for r in rows if not r['in_prior_GENCODE37_run']],indent=2))
