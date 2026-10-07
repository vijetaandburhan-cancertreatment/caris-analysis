"""Keep all completed public benchmark evidence in a verified, non-iCloud resident copy."""
from pathlib import Path
import hashlib,json,subprocess,shutil,datetime
P=Path(__file__).resolve().parent
D=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct5-fusion-robustness-v1')
assert (P/'execution-complete.json').exists()
assert (P/'results.json').exists()
assert json.loads((P/'results.json').read_text())['runs_completed']==45
assert shutil.disk_usage(D.parent).free>3*1024**3
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
stage=json.loads((P/'immutable-resident-stage.json').read_text())
assert stage['all_SHA256_match'] and stage['all_resident_flags_not_dataless']
verified_immutable={r['relative_path']:r for r in stage['files']}
review=json.loads((P/'independent-final-claims-review.json').read_text())
assert review['status']=='PASS' and review['runs_checked']==45
for filename,key in [('results.json','primary_results_sha256'),('summary.json','primary_summary_sha256'),('findings.txt','reviewed_findings_sha256')]:assert sha(P/filename)==review[key]
files=[f for f in P.rglob('*') if f.is_file() and not f.is_symlink() and '__pycache__' not in f.parts and 'audit-parser-history' not in f.parts and 'analysis-cache-initial-line-parser' not in f.parts and f.name not in ['resident-mirror-receipt.json']]
manifest={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_root':str(P),'resident_root':str(D),'scope':'Prespecified public-fixture benchmark, all inputs and outcomes, initial terminal-state halt, reviewed amendments, code and independent audits. No patient data or genome/index copies.','files':[]}
priority={'findings.txt':0,'summary.json':1,'summary.tsv':2,'results.json':3,'results.tsv':4,'independent-final-claims-review.json':5,'independent-outcome-audit.json':6}
for f in sorted(files,key=lambda f:(priority.get(str(f.relative_to(P)),100),str(f))):
    rel=f.relative_to(P);target=D/rel;target.parent.mkdir(parents=True,exist_ok=True)
    prior=verified_immutable.get(str(rel))
    want=prior['sha256'] if prior else sha(f)
    if prior:assert f.stat().st_size==prior['bytes']
    if target.exists():assert sha(target)==want,'Unexpected resident conflict '+str(target)
    else:
        cp=subprocess.run(['/bin/cp','-c','-p',str(f),str(target)],capture_output=True)
        if cp.returncode:shutil.copy2(f,target)
    assert sha(target)==want
    flags=subprocess.check_output(['stat','-f','%Sf',str(target)],text=True).strip();assert 'dataless' not in flags
    manifest['files'].append({'relative_path':str(rel),'bytes':f.stat().st_size,'sha256':want,'destination_flags':flags,'source_hash_provenance':'verified immutable-resident-stage.json' if prior else 'fresh source SHA256'})
manifest['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();manifest['file_count']=len(manifest['files']);manifest['total_logical_bytes']=sum(r['bytes'] for r in manifest['files'])
(D/'residency-manifest.json').write_text(json.dumps(manifest,indent=2))
receipt={'destination':str(D),'file_count':manifest['file_count'],'total_logical_bytes':manifest['total_logical_bytes'],'manifest_sha256':sha(D/'residency-manifest.json'),'all_destination_SHA256_match':True,'all_destination_flags_non_dataless':True,'finished_utc':manifest['finished_utc']}
(P/'resident-mirror-receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))
