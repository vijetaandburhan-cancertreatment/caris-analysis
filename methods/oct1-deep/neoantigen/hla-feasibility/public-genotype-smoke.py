"""Public upstream BAM control only; never runs patient genotype."""
from pathlib import Path
import json,subprocess,os,sys,time,resource
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
BASE=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/hla-tools')
ARC=BASE/'arcasHLA';PUB=BASE/'public-genotype-smoke';PUB.mkdir(exist_ok=True)
TMP=BASE/'temp';TMP.mkdir(exist_ok=True)
KAL=ROOT/'work/oct1-deep/fusion-next-step/public/kallisto-v0.44.0/kallisto'
env=os.environ.copy();env['PATH']=str(KAL)+os.pathsep+env['PATH']
args=[sys.executable,str(ARC/'scripts/genotype.py'),str(BASE/'smoke/public_test.R1.fq.gz'),str(BASE/'smoke/public_test.R2.fq.gz'),'-g','A,B,C,DPB1,DQB1,DQA1,DRB1','-o',str(PUB),'--temp',str(TMP),'-t','1','-v']
t=time.time()
with (OUT/'public-genotype-smoke.log').open('w') as log:
    p=subprocess.run(args,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
assert p.returncode==0,p.returncode
files=list(PUB.glob('*.genotype.json'));assert len(files)==1
got=json.loads(files[0].read_text());expected=json.loads((ARC/'test/expected_output/test.genotype.json').read_text())
def two(x):return sorted(':'.join(a.split(':')[:2]) for a in x)
comparison={k:{'observed':got.get(k),'official_expected':v,'two_field_match':k in got and two(v)==two(got[k])} for k,v in expected.items()}
status=all(z['two_field_match'] for z in comparison.values())
result={'status':'passed' if status else 'discordant_public_control','elapsed_seconds':time.time()-t,'comparison':comparison,'arcas_commit':'9fa54a212d134b0d9894d1fc19ec1bdc6f62eb55','IMGT_release':'3.46.0','kallisto_version':'0.44.0','note':'Only upstream public example processed. Documentation control uses older IMGT3.24; compare at two-field resolution explicitly. Public BAM has no quality strings; extractor supplied public-only Q40 placeholders for kallisto sequence-only operation. Patient extraction never uses placeholders.','full_command':args,'reference_directory_bytes':sum(p.stat().st_size for p in ARC.rglob('*') if p.is_file()),'source_python_peak_RSS_bytes_from_reference_build':1239826432}
final_expected=json.loads((ARC/'test/expected_output/test.partial_genotype.json').read_text())
result['final_typing_control_comparison']={k:{'observed':got.get(k),'official_final_expected':v,'two_field_match':k in got and two(v)==two(got[k])} for k,v in final_expected.items()}
result['final_typing_control_all_two_field_match']=all(z['two_field_match'] for z in result['final_typing_control_comparison'].values())
result['control_interpretation']='Reference version differs; standard-stage comparison and final partial-stage expectations both retained. Not a same-reference regression test, and any mismatch remains explicit.'
(OUT/'public-genotype-smoke.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
