import pathlib,sys,json,subprocess,gzip,hashlib,pysam,datetime,shutil
P=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(P));from align import run,STAR
C=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/genome-fusion-reference/compact_BCR_ABL1');OUT=P/'compact-controls';OUT.mkdir(exist_ok=True)
reads=[P/'public-controls'/f'arriba.R{x}.fastq.gz' for x in [1,2]];counts=[]
for f in reads:
 with gzip.open(f,'rt') as s:counts.append(sum(1 for _ in s)//4)
assert counts[0]==counts[1]
meta={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'control':'Official Arriba BCR::ABL1 synthetic FASTQs','control_pairs':counts[0],'reference':'Coordinate-preserving fullABL1/BCR +50kbflanks; other lociN/absent; not whole-genome equivalence','commands':[],'runs':{}}
for d in [1,8]:
 idx=C/f'index_D{d}_SA12';idx.mkdir(exist_ok=True);cmd=[str(STAR),'--runMode','genomeGenerate','--runThreadN','1','--genomeDir',str(idx),'--genomeFastaFiles',str(C/'assembly.fa'),'--sjdbGTFfile',str(C/'annotation.gtf'),'--sjdbOverhang','160','--genomeSAsparseD',str(d),'--genomeSAindexNbases','12','--limitGenomeGenerateRAM','1500000000','--outFileNamePrefix',str(idx/'build.')];meta['commands'].append(cmd)
 assert shutil.disk_usage(C).free>4*1024**3,'Need4GiB free before compact control'
 with (idx/'driver.log').open('w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
 meta['runs'][str(d)]=run(idx,C/'assembly.fa',C/'annotation.gtf',reads[0],reads[1],OUT/f'D{d}',threads=2,save_bam=True,expected_pairs=counts[0]);(OUT/'progress.json').write_text(json.dumps(meta,indent=2))

def load(d):
 out=OUT/f'D{d}';lines=(out/'fusions.tsv').read_text().splitlines();head=lines[0].lstrip('#').split('\t');fusions=[dict(zip(head,line.split('\t'))) for line in lines[1:]]
 with pysam.AlignmentFile(out/'Aligned.control.bam','rb') as bam:
  recs=sorted((r.query_name,r.flag,r.reference_name,r.reference_start,r.cigarstring,r.next_reference_name,r.next_reference_start,r.template_length,r.query_sequence,tuple(sorted(r.get_tags()))) for r in bam)
 junctions=sorted(l for l in (out/'STAR.Chimeric.out.junction').read_text().splitlines() if not l.startswith('#'))
 return fusions,recs,junctions
f1,r1,j1=load(1);f8,r8,j8=load(8)
check={'D1_positive':any({f['gene1'],f['gene2']}=={'BCR','ABL1'} and f['confidence']=='high' for f in f1),'D8_positive':any({f['gene1'],f['gene2']}=={'BCR','ABL1'} and f['confidence']=='high' for f in f8),'fusion_rows_identical':f1==f8,'all_alignment_tuples_identical':r1==r8,'chimeric_records_identical':j1==j8,'D1_alignment_records':len(r1),'D8_alignment_records':len(r8),'D1_chimeric_records':len(j1),'D8_chimeric_records':len(j8)}
meta.update(status='passed' if all(check[k] for k in ['D1_positive','D8_positive','fusion_rows_identical','all_alignment_tuples_identical','chimeric_records_identical']) else 'review_required',comparison=check,D1_fusions=f1,D8_fusions=f8,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());(OUT/'results.json').write_text(json.dumps(meta,indent=2));print(json.dumps(check,indent=2));assert meta['status']=='passed'
