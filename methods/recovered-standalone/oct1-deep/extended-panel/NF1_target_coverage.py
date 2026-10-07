"""GC/capture-uncalibrated descriptive coding depth, not copy-number calls."""
import pathlib,gzip,re,collections,json,csv,statistics,time,hashlib
import pysam
ROOT=pathlib.Path(__file__).resolve().parents[3];OUT=pathlib.Path(__file__).resolve().parent
SRC=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853');OLD=ROOT/'work/oct1-analysis'
GENES={'NF1'}
GTF=OUT.parent/'genomics/gencode.v37.annotation.gtf.gz'
def merge(a):
 out=[]
 for s,e in sorted(a):
  if out and s<=out[-1][1]:out[-1][1]=max(e,out[-1][1])
  else:out.append([s,e])
 return out
def write_tsv(name,rows):
 if not rows:return
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
genes={};tx={};lines=0
with gzip.open(GTF,'rt') as f:
 for ln in f:
  lines+=1
  if ln.startswith('#'):continue
  a=ln.rstrip().split('\t')
  if a[2] not in {'gene','transcript','CDS','exon'}:continue
  gene=re.search(r'gene_name "([^"]+)"',a[8]);gene=gene.group(1) if gene else None
  near9=False
  if gene not in GENES and not near9:continue
  if a[0] not in {'chr'+str(i) for i in range(1,23)}|{'chrX','chrY'}:continue
  attrs=collections.defaultdict(list)
  for k,val in re.findall(r'(\w+) "([^"]+)"',a[8]):attrs[k].append(val)
  if a[2]=='gene':
   genes[gene]={'gene':gene,'chrom':a[0],'gene_start0':int(a[3])-1,'gene_end0':int(a[4]),'strand':a[6],'gene_id':attrs['gene_id'][0],'is_9p21_neighborhood':near9}
  if 'transcript_id' not in attrs:continue
  tid=attrs['transcript_id'][0]
  t=tx.setdefault(tid,{'gene':gene,'chrom':a[0],'strand':a[6],'transcript_id':tid,'cds':[],'exons':[],'tags':[]})
  t['tags']=list(set(t['tags']+attrs['tag']))
  if a[2]=='CDS':t['cds'].append([int(a[3])-1,int(a[4])])
  if a[2]=='exon':t['exons'].append([int(a[3])-1,int(a[4])])
def priority(t):
 tags=t['tags'];ap=[int(z[-1]) for z in tags if z.startswith('appris_principal_')]
 return ('MANE_Select' in tags,-min(ap) if ap else -99,'basic' in tags,sum(e-s for s,e in merge(t['cds'])))
regions=[]
for gene,g in genes.items():
 ts=[t for t in tx.values() if t['gene']==gene and t['cds']]
 if not ts:continue
 selected=max(ts,key=priority);union=merge([z for t in ts for z in t['cds']])
 g.update(coding_transcript_count=len(ts),selected_transcript=selected['transcript_id'],selected_transcript_tags=selected['tags'],selected_CDS=merge(selected['cds']),all_transcript_CDS_union=union)
 for i,(s,e) in enumerate(union):regions.append((gene,g['chrom'],s,e,i))
(OUT/'NF1-target-annotation.json').write_text(json.dumps({'source':'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz','sha256':hashlib.sha256(GTF.read_bytes()).hexdigest(),'gtf_lines_streamed_to_EOF':lines,'coordinate_system':'0-based half-open','selection':'Prioritize MANE_Select, then APPRIS principal lowest rank, basic, then longest CDS. All-transcript union separately evaluated. Not Caris assay target BED.','genes':genes},indent=2)+'\n')

def callback(r):return not r.flag & (4|256|512|1024|2048) and r.mapping_quality>=20
def stats(d):
 return {'bases':len(d),'mean_read_depth':round(statistics.mean(d),3) if d else None,'median_read_depth':statistics.median(d) if d else None,'min_read_depth':min(d) if d else None,'max_read_depth':max(d) if d else None,**{f'bases_ge_{n}':sum(z>=n for z in d) for n in [1,10,20,50,100,200]}}
intervalrows=[];summary=[];lows=[];started=time.time()
with pysam.AlignmentFile(str(SRC/'DNA_TN26-279853.bam'),'rb',index_filename=str(OLD/'DNA_TN26-279853.bam.bai'),threads=1) as bam:
 for gene,g in genes.items():
  if 'selected_CDS' not in g:continue
  depth={}
  for i,(s,e) in enumerate(g['all_transcript_CDS_union']):
   cov=bam.count_coverage(g['chrom'],s,e,quality_threshold=20,read_callback=callback)
   d=[sum(z) for z in zip(*cov)];depth.update(zip(range(s,e),d))
   intervalrows.append({'gene':gene,'chrom':g['chrom'],'start0':s,'end0':e,'union_interval_index':i,**stats(d)})
   bad=merge([(s+j,s+j+1) for j,v in enumerate(d) if v<20])
   for a,b in bad:lows.append({'gene':gene,'chrom':g['chrom'],'start0':a,'end0':b,'length':b-a,'mean_read_depth':round(statistics.mean(depth[p] for p in range(a,b)),3)})
  for name,rs in [('all_transcript_CDS_union',g['all_transcript_CDS_union']),('selected_transcript_CDS',g['selected_CDS'])]:
   d=[depth[p] for s,e in rs for p in range(s,e)]
   summary.append({'gene':gene,'chrom':g['chrom'],'scope':name,'selected_transcript':g['selected_transcript'],**stats(d)})
  print(gene,summary[-1]['bases'],summary[-1]['median_read_depth'],summary[-1]['bases_ge_20'],round(time.time()-started,1),flush=True)
write_tsv('NF1-coding-coverage-summary.tsv',summary);write_tsv('NF1-coding-coverage-intervals.tsv',intervalrows);write_tsv('NF1-coding-coverage-low20-runs.tsv',lows)
(OUT/'NF1-coding-coverage-method.json').write_text(json.dumps({'started_unix':started,'elapsed_seconds':time.time()-started,'pysam':pysam.__version__,'method':'pysam count_coverage per merged GENCODE v37 coding intervals; A/C/G/T aligned bases; BQ>=20, MAPQ>=20, exclude unmapped/secondary/supplementary/QCfail/duplicate flags. Overlapping mates counted separately; read depth, not independent molecule depth. Deletion bases/Ns not counted. No full-gene or genome search for new variants.','limitations':['Public CDS is not the actual assay capture BED.','Coverage thresholds are descriptive, not validated clinical sensitivity or somatic-callability thresholds.','No calibrated normal/GC/capture normalization, purity/ploidy fit or segmented absolute copy numbers. Never infer gene deletion state or boundaries directly from these depths.','Alternative CDS union may include transcripts not intended by clinical assay.','BAP1 deletion itself causes reduced base coverage over the altered sequence; this is not an exon copy deletion.']},indent=2)+'\n')
