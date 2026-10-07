"""Disable population priors with the exact same local FASTQs.

Initial saved-alignment replay failed in the upstream six-vs-nine-field loader;
that error is preserved separately. Do not modify official inference source.
"""
from pathlib import Path
import json, os, subprocess, sys, time, signal
OUT=Path(__file__).resolve().parent
runs=json.loads((OUT/'patient-genotype-runs.json').read_text());assert len(runs)==2
def two(x):return sorted(':'.join(a.split(':')[:2]) for a in x)
result=[]
env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
kal=OUT.parents[3]/'work/oct1-deep/fusion-next-step/public/kallisto-v0.44.0/kallisto'
env['PATH']=str(kal)+os.pathsep+env['PATH']
for run in runs:
    original=run['command'];output=Path(run['output_dir']+'-no-prior-fastq');output.mkdir(exist_ok=False)
    args=original.copy()
    args[args.index('-p')+1]='none';args[args.index('-o')+1]=str(output)
    temp=output/'temp';temp.mkdir();args[args.index('--temp')+1]=str(temp)
    logfile=OUT/('patient-'+run['version']+'-no-prior-fastq.arcas.log');args[args.index('--log')+1]=str(logfile)
    start=time.time();peak=0;guard=None
    with (OUT/('patient-'+run['version']+'-no-prior-fastq.driver.log')).open('w') as f:
        p=subprocess.Popen(args,stdout=f,stderr=subprocess.STDOUT,env=env,start_new_session=True)
        while p.poll() is None:
            rows=[list(map(int,a.split())) for a in subprocess.check_output(['ps','-axo','pid=,ppid=,rss='],text=True).splitlines()]
            ids={p.pid};changed=True
            while changed:
                prev=len(ids);ids.update(pid for pid,ppid,rss in rows if ppid in ids);changed=len(ids)>prev
            rss=sum(rss*1024 for pid,ppid,rss in rows if pid in ids);peak=max(peak,rss)
            if rss>6*1024**3:guard='RSS cap'
            elif time.time()-start>600:guard='ten minute sensitivity cap'
            if guard:os.killpg(p.pid,signal.SIGTERM);p.wait(timeout=20);break
            time.sleep(1)
        rc=p.wait()
    d={'version':run['version'],'returncode':rc,'guard':guard,'elapsed_seconds':time.time()-start,'peak_process_tree_RSS_bytes_1sec_sampled':peak,'command':args}
    if rc==0 and not guard:
        observed=json.loads(next(output.glob('*.genotype.json')).read_text())
        d['genotypes_no_prior']=observed;d['comparison']={gene:{'default_prior':alleles,'no_prior':observed.get(gene),'two_field_match':gene in observed and two(alleles)==two(observed[gene])} for gene,alleles in run['genotypes'].items()}
    result.append(d);(OUT/'patient-prior-sensitivity.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(d),flush=True)
    assert rc==0 and not guard
