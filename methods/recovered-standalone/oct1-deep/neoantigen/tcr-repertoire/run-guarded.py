"""Task-local TRUST4 stages with resource/disk/elapsed guards and provenance."""
from pathlib import Path
import json,os,signal,subprocess,sys,time,shutil
OUT=Path(__file__).resolve().parent
BASE=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/trust4-tools')
label,cwd,*command=sys.argv[1:];assert label and command
before=sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())
started=time.time();peak=0;guard=None;checkpoints=[]
with (OUT/(label+'.log')).open('w') as log:
    p=subprocess.Popen(command,cwd=cwd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
    while p.poll() is None:
        rows=[tuple(map(int,line.split())) for line in subprocess.check_output(['ps','-axo','pid=,ppid=,rss='],text=True).splitlines()]
        ids={p.pid};changed=True
        while changed:
            old=len(ids);ids.update(pid for pid,ppid,rss in rows if ppid in ids);changed=len(ids)>old
        rss=sum(rss*1024 for pid,ppid,rss in rows if pid in ids);peak=max(peak,rss)
        size=sum(x.stat().st_size for x in BASE.rglob('*') if x.is_file())
        free=shutil.disk_usage(BASE).free;elapsed=time.time()-started
        if rss>4*1024**3:guard='4GiB process-tree RSS cap'
        elif size>1024**3:guard='1GiB total task additional disk cap'
        elif free<2*1024**3:guard='2GiB free disk guard'
        elif elapsed>1800:guard='30 minute stage guard'
        status={'stage':label,'elapsed_seconds':elapsed,'process_tree_RSS_bytes':rss,'process_tree_peak_RSS_bytes':peak,'directory_bytes':size,'free_disk_bytes':free,'guard':guard}
        (OUT/'progress.json').write_text(json.dumps(status,indent=2)+'\n')
        if len(checkpoints)==0 or elapsed-checkpoints[-1]['elapsed_seconds']>30:checkpoints.append(status)
        if guard:
            os.killpg(p.pid,signal.SIGTERM);p.wait(timeout=20);break
        time.sleep(2)
    rc=p.wait()
result={'label':label,'cwd':cwd,'command':command,'returncode':rc,'guard':guard,'elapsed_seconds':time.time()-started,'process_tree_peak_RSS_bytes_2sec_sampled':peak,'directory_bytes_before':before,'directory_bytes_after':sum(x.stat().st_size for x in BASE.rglob('*') if x.is_file()),'checkpoints':checkpoints,'status':'completed' if rc==0 and guard is None else 'failed_or_guarded_stop'}
(OUT/(label+'-result.json')).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2));sys.exit(0 if result['status']=='completed' else 1)
