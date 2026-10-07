"""Official public control on its documented IMGT3.24 version. No patient input."""
from pathlib import Path
import shutil,urllib.request,json,hashlib,sys,os,subprocess,time
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
BASE=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/hla-tools')
ARC=BASE/'arcasHLA';CTRL=BASE/'arcasHLA-control-3.24';REF=CTRL/'dat/IMGTHLA'
COMMIT='c5acf7a4342869351b2382b1cc1d1b5763e7e04e'
KAL=ROOT/'work/oct1-deep/fusion-next-step/public/kallisto-v0.44.0/kallisto'
os.environ['PATH']=str(KAL)+os.pathsep+os.environ['PATH']
CTRL.mkdir(exist_ok=True);shutil.copytree(ARC/'scripts',CTRL/'scripts',dirs_exist_ok=True);shutil.copytree(ARC/'dat/info',CTRL/'dat/info',dirs_exist_ok=True);(CTRL/'dat/ref').mkdir(exist_ok=True)
manifest=[]
for n in ['hla.dat','wmda/hla_nom_g.txt','wmda/hla_nom_p.txt']:
    p=REF/n;p.parent.mkdir(parents=True,exist_ok=True);url=f'https://raw.githubusercontent.com/ANHIG/IMGTHLA/{COMMIT}/{n}'
    if not p.exists():
        with urllib.request.urlopen(url,timeout=60) as r,p.open('wb') as f:shutil.copyfileobj(r,f)
    manifest.append({'path':str(p),'url':url,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
print('Downloaded same-version control resources',flush=True)
sys.path.insert(0,str(CTRL/'scripts'));import reference
reference.hla_dat_version=lambda print_version=False:COMMIT
ow=reference.write_reference
def main_only(seq,info,fasta,idx,database,typ):
    if typ=='partial':return
    return ow(seq,info,fasta,idx,database,typ)
reference.write_reference=main_only
t=time.time();reference.build_convert(False);reference.build_fasta()
print('Control reference built',flush=True)
PUB=BASE/'public-control-3.24';PUB.mkdir(exist_ok=True);TMP=BASE/'temp';TMP.mkdir(exist_ok=True)
args=[sys.executable,str(CTRL/'scripts/genotype.py'),str(BASE/'smoke/public_test.R1.fq.gz'),str(BASE/'smoke/public_test.R2.fq.gz'),'-g','A,B,C,DPB1,DQB1,DQA1,DRB1','-o',str(PUB),'--temp',str(TMP),'-t','1','-v']
with (OUT/'same-version-control-genotype.log').open('w') as f:p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,timeout=180)
assert p.returncode==0,p.returncode
got=json.loads(next(PUB.glob('*.genotype.json')).read_text());expected=json.loads((ARC/'test/expected_output/test.genotype.json').read_text())
def two(v):return sorted(':'.join(a.split(':')[:2]) for a in v)
cmp={k:{'observed':got.get(k),'official_expected':v,'two_field_match':k in got and two(got[k])==two(v)} for k,v in expected.items()}
result={'status':'passed' if all(z['two_field_match'] for z in cmp.values()) else 'failed','IMGT_version':'3.24.0','elapsed_seconds':time.time()-t,'comparison':cmp,'public_reference_manifest':manifest,'total_public_tools_reference_directory_bytes':sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file()),'note':'Only documented upstream public dataset. Standard complete-allele stage compared with standard-stage official expected values at same reference version. No patient genotype.'}
(OUT/'same-version-control.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
