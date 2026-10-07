"""Single-sample descriptive fragment density; NOT a calibrated CN caller.

Protein-coding GENCODE37 exon masks are an explicit proxy for capture space.
Aligned proper-pair read1 midpoints are counted once per query name by flags;
duplicate-flagged, secondary/supplementary/QC-fail records are excluded.
"""
import pathlib,sys,gzip,re,collections,json,time,bisect,hashlib,resource
import numpy as np,pysam,py2bit
OUT=pathlib.Path(__file__).resolve().parent; ROOT=OUT.parents[2]
SRC=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
GTF=ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz'
EXCLUDE=4|8|256|512|1024|2048
def merge(a):
 out=[]
 for s,e in sorted(a):
  if out and s<=out[-1][1]:out[-1][1]=max(out[-1][1],e)
  else:out.append([s,e])
 return out
started=time.time();exons=collections.defaultdict(list);genes={};cds=collections.defaultdict(list)
with gzip.open(GTF,'rt') as f:
 for ln in f:
  if ln.startswith('#'):continue
  a=ln.rstrip().split('\t')
  if a[0] not in {'chr'+str(i) for i in range(1,23)}|{'chrX','chrY'}:continue
  if 'gene_type "protein_coding"' not in a[8]:continue
  if a[2] not in {'gene','exon','CDS'}:continue
  gene=re.search('gene_name "([^"]+)"',a[8]).group(1)
  if a[2]=='gene':genes[gene]={'chrom':a[0],'start0':int(a[3])-1,'end0':int(a[4]),'strand':a[6]}
  if a[2]=='exon':exons[a[0]].append([int(a[3])-1,int(a[4])])
  if a[2]=='CDS':cds[gene].append([int(a[3])-1,int(a[4])])
for ch in exons:exons[ch]=merge(exons[ch])
for gene in cds:cds[gene]=merge(cds[gene])
(OUT/'annotation.json').write_text(json.dumps({'coordinate_system':'0-based half-open','GENCODE_version':37,'genes':genes,'protein_coding_exon_union':exons,'protein_coding_CDS_union_by_gene':cds},separators=(',',':'))+'\n')
tb=py2bit.open(str(OUT/'public/hg38.2bit'),storeMasked=True)
bam=pysam.AlignmentFile(str(SRC/'DNA_TN26-279853.bam'),'rb',index_filename=str(ROOT/'work/oct1-analysis/DNA_TN26-279853.bam.bai'),threads=1)
summary={};(OUT/'genome-bins').mkdir(exist_ok=True)
for chrom in ['chr'+str(i) for i in range(1,23)]+['chrX','chrY']:
 t=time.time();L=bam.get_reference_length(chrom);N=(L+999)//1000
 dest=OUT/'genome-bins'/f'{chrom}.npz'
 if dest.exists():print(chrom,'already complete',flush=True);continue
 # 1-byte static mask: bit0 exon+100, bit1 exon+500, bit2 exon+2000.
 mask=np.zeros(L,dtype=np.uint8)
 for pad,bit in [(100,1),(500,2),(2000,4)]:
  for s,e in exons[chrom]:mask[max(0,s-pad):min(L,e+pad)] |= bit
 seq=tb.sequence(chrom).encode('ascii');a=np.frombuffer(seq,dtype=np.uint8)
 start=np.arange(0,L,1000);gc=((a==ord('G'))|(a==ord('C'))|(a==ord('g'))|(a==ord('c')))
 acgt=(a==65)|(a==67)|(a==71)|(a==84)|(a==97)|(a==99)|(a==103)|(a==116)
 soft=(a>=97)&(a<=122)
 arrays={'length':np.minimum(1000,L-start).astype(np.int32),'acgt_bases':np.add.reduceat(acgt.astype(np.int32),start),'gc_bases':np.add.reduceat(gc.astype(np.int32),start),'repeat_bases':np.add.reduceat(soft.astype(np.int32),start)}
 # Denominators and GC restricted to the same off-exon spaces as fragment centers.
 for pad,bit in [(100,1),(500,2),(2000,4)]:
  keep=(mask&bit)==0
  arrays[f'off{pad}_acgt_bases']=np.add.reduceat((keep&acgt).astype(np.int32),start)
  arrays[f'off{pad}_gc_bases']=np.add.reduceat((keep&gc).astype(np.int32),start)
  arrays[f'off{pad}_repeat_bases']=np.add.reduceat((keep&soft).astype(np.int32),start)
 for key in ['all','off100','off500','off2000','off500_mq60']:
  arrays[f'fragment_midpoints_{key}']=np.zeros(N,dtype=np.int32)
 del seq,a,gc,acgt,soft,keep
 ctr=collections.Counter();tlens=collections.Counter()
 for r in bam.fetch(chrom):
  ctr['alignment_records_seen']+=1
  if not r.is_read1 or r.flag&EXCLUDE or not r.is_proper_pair:continue
  if r.mapping_quality<30 or r.next_reference_id!=r.reference_id:continue
  if not 50<=abs(r.template_length)<=1000:continue
  left=min(r.reference_start,r.next_reference_start);mid=left+abs(r.template_length)//2
  if not 0<=mid<L:continue
  ctr['accepted_read1_fragments']+=1;tlens[abs(r.template_length)]+=1;i=mid//1000;m=mask[mid]
  arrays['fragment_midpoints_all'][i]+=1
  if not m&1:arrays['fragment_midpoints_off100'][i]+=1
  if not m&2:
   arrays['fragment_midpoints_off500'][i]+=1
   if r.mapping_quality>=60:arrays['fragment_midpoints_off500_mq60'][i]+=1
  if not m&4:arrays['fragment_midpoints_off2000'][i]+=1
 np.savez_compressed(dest,**arrays)
 result={'counts':dict(ctr),'template_lengths':dict(tlens),'elapsed_seconds':time.time()-t,'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 (OUT/'genome-bins'/f'{chrom}.json').write_text(json.dumps(result,indent=2)+'\n');summary[chrom]=result
 print(chrom,dict(ctr),round(time.time()-t,1),'RSS',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,flush=True)
 del mask,arrays
(OUT/'genome-fragment-profile.method.json').write_text(json.dumps({'elapsed_seconds':time.time()-started,'pysam':pysam.__version__,'py2bit':py2bit.__version__,'bins_bp':1000,'fragment_filters':{'read1':True,'proper_pair':True,'exclude_flag_mask':EXCLUDE,'MAPQ_min':30,'insert_size_min':50,'insert_size_max':1000,'mate_MAPQ':'not independently filtered; MAPQ60 sensitivity on read1 only','PCR':'duplicate flags excluded; no UMI molecule assertion'},'capture_proxy':'GENCODE37 all protein-coding transcript exon union with100/500/2000bp flanks, NOT actual baits','limitations':['Unnormalized single-sample fragment density is NOT copy number.','Midpoints proxy physical fragments; cannot establish unique original molecules.','GC and repeat adjustment cannot remove unknown bait/capture/protocol biases.','Reference is UCSC hg38 primary rather than exact private masked Caris FASTA; unplaced/alternate contigs excluded.','Absolute copy number/purity/ploidy and somatic LOH are not identifiable without additional assumptions or independent reference.']},indent=2)+'\n')
