"""Local research prediction. Empty terminal flanks are retained as empty strings."""
import os,csv,json,hashlib,time,importlib.metadata
from pathlib import Path
os.environ['MHCFLURRY_DATA_DIR']='/Users/burhanazeem/.local/share/codex/caris-analysis/mhcflurry-models'
os.environ['MHCFLURRY_DEVICE']='cpu'
os.environ['CUDA_VISIBLE_DEVICES']=''
import torch
torch.set_num_threads(4)
from mhcflurry import Class1PresentationPredictor
OUT=Path(__file__).resolve().parent
rows=list(csv.DictReader((OUT/'mhcflurry-input.csv').open()))
assert all(isinstance(r[k],str) for r in rows for k in ['peptide','n_flank','c_flank'])
sets={r['allele']:r['allele'].split(';') for r in rows}
p=Class1PresentationPredictor.load()
print('Loaded',importlib.metadata.version('mhcflurry'),'rows',len(rows),flush=True)
start=time.monotonic()
pred=p.predict(peptides=[r['peptide'] for r in rows],alleles=sets,sample_names=[r['allele'] for r in rows],n_flanks=[r['n_flank'] for r in rows],c_flanks=[r['c_flank'] for r in rows],include_affinity_percentile=True,verbose=0)
assert len(pred)==len(rows)
for name in ['candidate_id','variant','gene','sequence_type','allele','normal_reference_exact_match_count']:
    pred[name]=pred['peptide_num'].map({i:r[name] for i,r in enumerate(rows)})
for _,x in pred.iterrows():assert rows[x['peptide_num']]['peptide']==x['peptide'] and rows[x['peptide_num']]['allele']==x['sample_name']
assert not pred['affinity'].isna().any()
pred.to_csv(OUT/'mhcflurry-predictions.csv',index=False)
modelroot=Path(os.environ['MHCFLURRY_DATA_DIR'])/'2.3.0/models_class1_presentation'
models=[{'path':str(f.relative_to(modelroot)),'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sorted(modelroot.rglob('*')) if f.is_file()]
meta={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'mhcflurry_version':importlib.metadata.version('mhcflurry'),'torch_version':torch.__version__,'model_release':'2.3.0','model_bundle_url':'https://github.com/openvax/mhcflurry/releases/download/2.3.0/models_class1_presentation.20260928.tar.bz2','official_archive_sha256':'768e360832b3d45dfe9d5e8cde4ffcb90abdb87773ca65f3a86ba7cc20f614db','elapsed_seconds':time.monotonic()-start,'output_rows':len(pred),'model_files':models,'fix':'Initial CLI read empty terminal flanks as NaN; rerun through documented Python API using csv module preserves true empty strings. No artificial residues added.','limits':['Research model scores are not observed binding, peptide display, T-cell response, safety or clinical efficacy.','Presentation percentile/score are model outputs, not probability of benefit.','A low IC50/rank may occur for germline or self peptides. Tumor specificity unestablished without matched normal.','Single mutation on specified reference isoform; actual full-length transcript and linked variants unmeasured.']}
(OUT/'mhcflurry-run-manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps({k:v for k,v in meta.items() if k!='model_files'},indent=2),flush=True)
# Keep support provenance and independently optimized allele summaries explicit.
import subprocess,sys
subprocess.run([sys.executable,str(OUT/'generate-candidate-summary.py')],check=True)
