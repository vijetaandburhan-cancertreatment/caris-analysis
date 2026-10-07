"""Two pinned references, same local patient input, sequential bounded jobs."""
from pathlib import Path
import json, os, signal, subprocess, sys, time, shutil
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
BASE=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/hla-tools')
KAL=ROOT/'work/oct1-deep/fusion-next-step/public/kallisto-v0.44.0/kallisto'
env=os.environ.copy();env['PATH']=str(KAL)+os.pathsep+env['PATH'];env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
extract=json.loads((OUT/'patient-extraction.json').read_text());assert extract['status']=='passed'
reads=[x['path'] for x in extract['files']]
GENES='A,B,C,DPA1,DPB1,DQA1,DQB1,DRB1,DRB5'
all_runs=[]
for version,source in [('3.24.0','arcasHLA-control-3.24'),('3.46.0','arcasHLA')]:
    dest=BASE/('patient-genotype-'+version);dest.mkdir(exist_ok=False)
    temp=BASE/('patient-temp-'+version);temp.mkdir(exist_ok=False)
    args=[sys.executable,str(BASE/source/'scripts/genotype.py'),*reads,'-g',GENES,'-o',str(dest),'--temp',str(temp),'-t','1','-p','prior','--min_count','75','--tolerance','0.000001','--max_iterations','1000','--drop_iterations','20','--drop_threshold','0.1','--zygosity_threshold','0.15','--log',str(OUT/('patient-'+version+'.arcas.log')),'-v']
    started=time.time();peak=0;guard=None
    print(json.dumps({'starting_version':version,'command':args}),flush=True)
    with (OUT/('patient-'+version+'.driver.log')).open('w') as log:
        p=subprocess.Popen(args,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        while p.poll() is None:
            rows=[]
            for line in subprocess.check_output(['ps','-axo','pid=,ppid=,rss='],text=True).splitlines():
                a=line.split()
                if len(a)==3: rows.append(tuple(map(int,a)))
            ids={p.pid};changed=True
            while changed:
                prior=len(ids);ids.update(pid for pid,ppid,rss in rows if ppid in ids);changed=len(ids)>prior
            rss=sum(rss*1024 for pid,ppid,rss in rows if pid in ids);peak=max(peak,rss)
            elapsed=time.time()-started;free=shutil.disk_usage(BASE).free
            if rss>6*1024**3:guard='6GiB process-tree RSS cap'
            elif free<2*1024**3:guard='2GiB free disk guard'
            elif elapsed>1200:guard='20 minute per-reference limit'
            (OUT/'patient-genotype-progress.json').write_text(json.dumps({'version':version,'elapsed_seconds':elapsed,'process_tree_peak_RSS_bytes':peak,'free_disk_bytes':free,'guard':guard},indent=2)+'\n')
            if guard:
                os.killpg(p.pid,signal.SIGTERM);p.wait(timeout=20);break
            time.sleep(2)
        rc=p.wait()
    result={'version':version,'returncode':rc,'guard':guard,'elapsed_seconds':time.time()-started,'process_tree_peak_RSS_bytes_2sec_sampled':peak,'command':args,'output_dir':str(dest),'files':[str(x) for x in dest.rglob('*') if x.is_file()],'arcas_commit':'9fa54a212d134b0d9894d1fc19ec1bdc6f62eb55','kallisto_version':'0.44.0','status':'completed' if rc==0 and not guard else 'failed_or_bounded_stop'}
    if result['status']=='completed':
        files=list(dest.glob('*.genotype.json'));assert len(files)==1
        result['genotypes']=json.loads(files[0].read_text())
    all_runs.append(result);(OUT/'patient-genotype-runs.json').write_text(json.dumps(all_runs,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    assert result['status']=='completed', result
