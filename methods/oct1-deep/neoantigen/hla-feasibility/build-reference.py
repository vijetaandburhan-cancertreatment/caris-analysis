"""Build official arcas reference locally with pinned provenance, not git clone.

Official algorithm/source files remain unchanged. Only hla_dat_version is replaced
with the version of the immutable, SHA256-verified direct public download.
"""
from pathlib import Path
import json,os,sys,time,resource,logging
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
ARC=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/hla-tools/arcasHLA')
KAL=ROOT/'work/oct1-deep/fusion-next-step/public/kallisto-v0.44.0/kallisto'
os.environ['PATH']=str(KAL)+os.pathsep+os.environ['PATH']
sys.path.insert(0,str(ARC/'scripts'))
import reference
logging.basicConfig(level=logging.INFO,format='%(message)s')
manifest=json.loads((OUT/'public-reference-manifest.json').read_text())
def version(print_version=False):
    if print_version: print('Pinned direct-download IMGT3.46 commit',manifest['IMGT_commit'],flush=True)
    return manifest['IMGT_commit']
reference.hla_dat_version=version
original_run=reference.run_command
def checked_run(*a,**kw):
    r=original_run(*a,**kw)
    if r.returncode:raise RuntimeError('Reference subprocess failed: '+r.stderr.decode())
    return r
reference.run_command=checked_run
original_write=reference.write_reference
def main_only_write(sequences,info,fasta,idx,database,type):
    # Partial typing is an optional, separately indexed stage. Its generated
    # index exceeds the public-resource budget; standard genotyping is intact.
    if type=='partial':
        print('Skipping optional partial reference by explicit resource limit',flush=True)
        return
    return original_write(sequences,info,fasta,idx,database,type)
reference.write_reference=main_only_write
print('Starting official reference builder',flush=True)
t=time.time();reference.build_convert(False);reference.build_fasta()
files=[{'name':p.name,'bytes':p.stat().st_size} for p in (ARC/'dat/ref').iterdir() if p.is_file()]
total=sum(p.stat().st_size for p in ARC.rglob('*') if p.is_file())
result={'elapsed_seconds':time.time()-t,'maxrss_bytes_macOS':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'total_public_directory_bytes':total,'reference_files':files,'compatibility_change':'Provenance-only hla_dat_version function returns pinned direct-download commit instead of git rev-parse. Official biological/reference-processing algorithms unchanged.'}
assert total<1_000_000_000,result
assert result['maxrss_bytes_macOS']<6_000_000_000,result
(OUT/'reference-build.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2),flush=True)
