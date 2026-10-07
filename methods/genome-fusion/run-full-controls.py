"""Public synthetic controls against full GENCODE37 primary assembly; no filter tuning."""
import pathlib,sys,json,gzip,datetime,hashlib
P=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(P))
from align import run
C=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/genome-fusion-reference')
OUT=P/'full-reference-controls';OUT.mkdir(exist_ok=True)
meta={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'index':'Full primary GENCODE37 D8/SAindex12/overhang160','runs':{},'controls':{}}
for name in ['arriba','starfusion']:
 reads=[P/'public-controls'/f'{name}.R{x}.fastq.gz' for x in [1,2]]
 counts=[]
 for f in reads:
  with gzip.open(f,'rt') as s: counts.append(sum(1 for _ in s)//4)
 assert counts[0]==counts[1]
 meta['runs'][name]=run(C/'STAR_G37_primary_D8_SA12',C/'GRCh38.primary_assembly.genome.fa',C/'gencode.v37.primary_assembly.annotation.gtf',*reads,OUT/name,threads=4,save_bam=True,expected_pairs=counts[0])
 lines=(OUT/name/'fusions.tsv').read_text().splitlines();head=lines[0].lstrip('#').split('\t')
 calls=[dict(zip(head,s.split('\t'))) for s in lines[1:]]
 meta['controls'][name]={'input_pairs':counts[0],'calls':calls,'accepted_gene_pairs':[[f['gene1'],f['gene2']] for f in calls]}
 if name=='arriba':
  positive=any({f['gene1'],f['gene2']}=={'BCR','ABL1'} and f['confidence']=='high' for f in calls)
  meta['BCR_ABL1_execution_control_pass']=positive
  assert positive,'Full-reference positive execution control failed'
 (OUT/'progress.json').write_text(json.dumps(meta,indent=2))
meta.update(status='complete',finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),limitations=['Known BCR::ABL1 synthetic control with known-fusion prior establishes execution, not general novel fusion/FFPE sensitivity.','The second published STAR-Fusion fixture is independently generated; differences between callers and annotations must be reviewed rather than forced by threshold tuning.','D8 compact controls lost one support read relative to D1; sparse-index sensitivity is not equivalent to standard index.'])
(OUT/'results.json').write_text(json.dumps(meta,indent=2));print(json.dumps({k:meta[k] for k in ['status','BCR_ABL1_execution_control_pass','finished_utc']},indent=2))
