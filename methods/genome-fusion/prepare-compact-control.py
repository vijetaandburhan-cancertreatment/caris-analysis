import pathlib,json,pysam,hashlib,datetime
P=pathlib.Path(__file__).resolve().parent;cache=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis');ref=cache/'genome-fusion-reference';out=ref/'compact_BCR_ABL1';out.mkdir(exist_ok=True)
regions={'chr9':('ABL1',130713016,130887675),'chr22':('BCR',23179704,23318037)};fa=ref/'GRCh38.primary_assembly.genome.fa'
if not pathlib.Path(str(fa)+'.fai').exists():pysam.faidx(str(fa))
fasta=pysam.FastaFile(str(fa));details=[]
with (out/'assembly.fa').open('w') as f:
 for chrom,(gene,start,end) in regions.items():
  lo=start-1-50000;hi=end+50000;f.write('>'+chrom+'\n')
  # Preserve actual GRCh38 coordinates; every base before retained flank is N.
  full,rem=divmod(lo,100)
  block=('N'*100+'\n')*10000
  for j in range(0,full,10000):f.write(block[:min(10000,full-j)*101])
  retained=fasta.fetch(chrom,lo,hi);seq='N'*rem+retained
  for j in range(0,len(seq),100):f.write(seq[j:j+100]+'\n')
  details.append({'chromosome':chrom,'gene':gene,'retained_start1':lo+1,'retained_end1':hi,'gene_start1':start,'gene_end1':end,'retained_sequence_sha256':hashlib.sha256(retained.encode()).hexdigest()})
count=0
with (ref/'gencode.v37.primary_assembly.annotation.gtf').open() as f,(out/'annotation.gtf').open('w') as g:
 for line in f:
  if not any(f'gene_name "{gene}";' in line for gene in ['ABL1','BCR']):continue
  x=line.split('\t');chrom=x[0];gene,a,b=regions[chrom];lo=a-50000;hi=b+50000
  assert lo<=int(x[3])<=int(x[4])<=hi;g.write(line);count+=1
meta={'generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'regions':details,'annotation_records':count,'only_full_ABL1_BCR_gene_annotation_retained':True,'coordinates':'Unshifted GRCh38 coordinates; longNprefix preserves coordinates. Other genomic loci absent or replaced byN. Not a specificity or whole-genome sensitivity benchmark.','files':{f.name:{'bytes':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in [out/'assembly.fa',out/'annotation.gtf']}}
(out/'manifest.json').write_text(json.dumps(meta,indent=2));(P/'compact-reference.json').write_text(json.dumps(meta,indent=2));print(json.dumps(meta,indent=2))
