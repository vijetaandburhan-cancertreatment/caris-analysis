"""Gate counter-only Arriba on frozen public BAMs, then re-run unchanged full RNA D8 mapping."""
import pathlib,sys,json,subprocess,hashlib,datetime,csv,shutil
P=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(P));import align
D=P/'diagnostics';BIN=D/'arriba-counter-build/arriba';CACHE=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis');R=CACHE/'genome-fusion-reference';RAW=CACHE/'TN26-279853'
def sha(f):
 with pathlib.Path(f).open('rb') as h:return hashlib.file_digest(h,'sha256').hexdigest()
def rows(f):
 ls=pathlib.Path(f).read_text().splitlines();rs=list(csv.DictReader([ls[0].lstrip('#')]+ls[1:],delimiter='\t'))
 for r in rs:
  if r.get('read_identifiers') not in ['.',None]:r['read_identifiers']=','.join(sorted(r['read_identifiers'].split(',')))
 return sorted(json.dumps(r,sort_keys=True) for r in rs)
meta={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'instrumented_binary_sha256':sha(BIN),'scope':'Diagnostic counters only, unchanged stock decision logic. No claim of fixing sparse-index sensitivity.','public_controls':{},'stock_outputs_frozen_sha256':{f:sha(P/'patient-D8'/f) for f in ['fusions.tsv','fusions.discarded.tsv','STAR.Chimeric.out.junction','STAR.SJ.out.tab','target-genes.bam']}}
(D/'rerun-manifest.json').write_text(json.dumps(meta,indent=2))
for name in ['arriba','starfusion']:
 original=P/'full-reference-controls'/name;out=D/'public-controls'/name;out.mkdir(parents=True,exist_ok=True)
 assert not (out/'fusions.tsv').exists(),'Refuse control overwrite'
 cmd=json.loads((original/'run.json').read_text())['Arriba_command'];cmd[0]=str(BIN)
 for arg,replacement in [('-x',str(original/'Aligned.control.bam')),('-o',str(out/'fusions.tsv')),('-O',str(out/'fusions.discarded.tsv'))]:cmd[cmd.index(arg)+1]=replacement
 with (out/'Arriba.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
 checks={f:{'byte_identical':sha(original/f)==sha(out/f),'all_rows_identical_after_row_and_readname_sort':rows(original/f)==rows(out/f),'stock_sha256':sha(original/f),'diagnostic_sha256':sha(out/f)} for f in ['fusions.tsv','fusions.discarded.tsv']}
 meta['public_controls'][name]={'command':cmd,'checks':checks};(D/'rerun-manifest.json').write_text(json.dumps(meta,indent=2))
 assert all(c['all_rows_identical_after_row_and_readname_sort'] for c in checks.values()),'Counter build changed public control output'
meta['public_controls_passed']=True;meta['patient_started_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(D/'rerun-manifest.json').write_text(json.dumps(meta,indent=2))
assert shutil.disk_usage(D).free>4*1024**3
align.ARRIBA=BIN
result=align.run(R/'STAR_G37_primary_D8_SA12',R/'GRCh38.primary_assembly.genome.fa',R/'gencode.v37.primary_assembly.annotation.gtf',RAW/'RNA_TN26-279853_S25.R1.fastq.gz',RAW/'RNA_TN26-279853_S25.R2.fastq.gz',D/'patient-D8-counter-rerun',threads=6,save_bam=False,expected_pairs=23209264)
meta['result']={'status':result['status'],'finished_utc':result['finished_utc'],'STAR_final':result['STAR_final']};(D/'rerun-manifest.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta['result'],indent=2))
