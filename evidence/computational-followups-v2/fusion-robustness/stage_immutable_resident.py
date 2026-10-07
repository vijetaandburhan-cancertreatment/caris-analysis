"""Parallel, bounded residency protection for CLOSED benchmark runs and frozen inputs only."""
from pathlib import Path
import hashlib,subprocess,json,shutil,concurrent.futures,datetime
P=Path(__file__).resolve().parent
D=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct5-fusion-robustness-v1')
assert json.loads((P/'execution-complete.json').read_text())['completed_runs']==45
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def copy(f):
    assert shutil.disk_usage(D.parent).free>3*1024**3
    relative=f.relative_to(P);t=D/relative;t.parent.mkdir(parents=True,exist_ok=True);want=sha(f)
    if t.exists():assert sha(t)==want
    else:
        result=subprocess.run(['/bin/cp','-c','-p',str(f),str(t)],capture_output=True)
        if result.returncode:shutil.copy2(f,t)
    assert sha(t)==want
    flags=subprocess.check_output(['stat','-f','%Sf',str(t)],text=True).strip();assert 'dataless' not in flags
    return {'relative_path':str(relative),'bytes':f.stat().st_size,'sha256':want,'resident_flags':flags}
files=[f for folder in ['runs','inputs'] for f in (P/folder).rglob('*') if f.is_file() and not f.is_symlink()]
rows=[];D.mkdir(parents=True,exist_ok=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=8) as e:
    futures=[e.submit(copy,f) for f in files]
    for result in concurrent.futures.as_completed(futures):
        rows.append(result.result())
        if len(rows)%50==0:print(f'Immutable resident copies verified: {len(rows)}/{len(files)}',flush=True)
out={'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'destination':str(D),'file_count':len(rows),'total_logical_bytes':sum(r['bytes'] for r in rows),'all_SHA256_match':True,'all_resident_flags_not_dataless':True,'max_workers':8,'files':sorted(rows,key=lambda r:r['relative_path'])}
(P/'immutable-resident-stage.json').write_text(json.dumps(out,indent=2));print(json.dumps({k:v for k,v in out.items() if k!='files'},indent=2))
