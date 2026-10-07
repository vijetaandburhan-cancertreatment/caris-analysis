"""Independent public common-SNP panel and query-name allele counts.

No genotypes or sequence data leave this machine. Public catalog selection is
predefined before reading patient allele counts; covered heterozygous markers
remain an ascertainment-limited tumor-only inference, not germline calls.
"""
import pathlib,json,bisect,csv,collections,hashlib,sys,time,statistics,gzip,resource
import numpy as np,pysam,py2bit
OUT=pathlib.Path(__file__).resolve().parent;ROOT=OUT.parents[2]
SRC=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
CHROMS=['chr3','chr5','chr9','chr17','chr2','chr8']
def write_tsv(path,rows):
 if not rows:return
 with path.open('w') as f:
  w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def merge(a):
 out=[]
 for s,e in sorted(a):
  if out and s<=out[-1][1]:out[-1][1]=max(e,out[-1][1])
  else:out.append([s,e])
 return out
def in_intervals(p,starts,ends):
 i=bisect.bisect_right(starts,p)-1
 return i>=0 and p<ends[i]
def select():
 ann=json.loads((OUT/'annotation.json').read_text());tb=py2bit.open(str(OUT/'public/hg38.2bit'))
 v=pysam.VariantFile(str(OUT/'public/1000G.GRCh38.20190312.sites.vcf.gz'))
 rows=[];summary={}
 for ch in CHROMS:
  regions=merge([(max(0,s-100),e+100) for s,e in ann['protein_coding_exon_union'][ch]])
  starts=[s for s,e in regions];ends=[e for s,e in regions]
  focus=[]
  for g in ['BAP1','RASA1','MTAP','CDKN2A','CDKN2B']:
   a=ann['genes'][g]
   if a['chrom']==ch:focus.append((max(0,a['start0']-250000),a['end0']+250000))
  focus=merge(focus);fs=[s for s,e in focus];fe=[e for s,e in focus]
  counts=collections.Counter();last_grid=-100000
  for r in v.fetch(ch.removeprefix('chr')):
   counts['public_records_seen']+=1
   if len(r.ref)!=1 or len(r.alts or [])!=1 or len(r.alts[0])!=1 or r.ref not in 'ACGT' or r.alts[0] not in 'ACGT':continue
   af=float(r.info['AF'][0]);sas=float(r.info['SAS_AF'][0]);p=r.pos-1
   if not .05<=af<=.95 or 'PASS' not in r.filter:continue
   exon=in_intervals(p,starts,ends);focal=in_intervals(p,fs,fe)
   grid=(.2<=af<=.8 and p-last_grid>=10000)
   if grid:last_grid=p
   if not(exon or focal or grid):continue
   ref=tb.sequence(ch,p,p+1).upper()
   if ref!=r.ref:counts['reference_mismatch_excluded']+=1;continue
   rows.append({'chrom':ch,'pos1':r.pos,'id':r.id or '.','ref':r.ref,'alt':r.alts[0],'public_AF':af,'public_SAS_AF':sas,'scope':'+'.join(k for k,b in [('exon100',exon),('focal250kb',focal),('grid10kb',grid)] if b)})
   counts['selected']+=1
  summary[ch]=dict(counts);print('panel',ch,dict(counts),flush=True)
 rows.sort(key=lambda r:(CHROMS.index(r['chrom']),r['pos1']))
 write_tsv(OUT/'public-common-SNP-panel.tsv',rows)
 for ch in CHROMS:
  with (OUT/f'public-common-SNP-panel.{ch}.bed').open('w') as f:
   for r in rows:
    if r['chrom']==ch:f.write(f"{ch}\t{r['pos1']-1}\t{r['pos1']}\n")
 (OUT/'public-panel-selection.json').write_text(json.dumps({'counts':summary,'reference':'1000G20190312 directly called GRCh38 biallelic panel; UCSC hg38 REF matches required','AF':'global alternate frequency0.05–0.95; SAS frequency recorded, not used to select','scope':'all protein-coding GENCODE37 exons+100bp on chr2/3/5/8/9/17, plus BAP1/RASA1/MTAP/CDKN2A/B±250kb all common SNPs, plus globalAF0.2–0.8 grid at least10kb apart; no patient allele evidence in selection','limitations':'These are population polymorphisms, not proof of heterozygosity or germline status in this individual. Population linkage and target enrichment violate independent-site assumptions.'},indent=2)+'\n')
 return rows
def command(bam,bed,out,index=None,region=None):
 args=['-B','-q','30','-Q','25','-d','100000','--ff','0xF0C','--rf','2','-x','--output-QNAME','--output-BP-5','--output-MQ','--no-output-ins','--no-output-ins','--no-output-del','--no-output-del','--no-output-ends','-l',str(bed),'-o',str(out)]
 if region:args+=['-r',region]
 if index:args+=['-X',str(bam),str(index)]
 else:args+=[str(bam)]
 return args
def parse_line(ln,marker):
 a=ln.rstrip('\n').split('\t')
 # With no FASTA, aligned bases are explicit. Output order from samtools is
 # bases/BQ, MQ, QNAME, BP5 (qualified below using synthetic control).
 if len(a)!=9:raise ValueError(f'Expected9fields,got{len(a)}')
 if int(a[3])==0:
  bases='';quals='';mq='';names=[];bp=[]
 else:
  bases=a[4];quals=a[5];mq=a[6];names=a[7].split(',');bp=[int(x) for x in a[8].split(',')]
 if any(len(x)!=len(bases) for x in [quals,mq,names,bp]):raise ValueError('Pileupparallelfieldlengthmismatch')
 groups=collections.defaultdict(list)
 for b,q,m,n,p in zip(bases,quals,mq,names,bp):
  b=marker['ref'] if b in '.,' else b.upper()
  groups[n].append((b,ord(q)-33,ord(m)-33,p))
 ctr=collections.Counter();read=collections.Counter();edge=collections.Counter();strand=collections.Counter()
 for b in bases:
  if b in 'ACGTacgt':strand[(b.upper(),'rev' if b.islower() else 'fwd')]+=1
 for n,items in groups.items():
  alleles={b for b,q,m,p in items if b in 'ACGT'}
  for b,q,m,p in items:read[b]+=1
  if len(alleles)>1:ctr['discordant']+=1;continue
  if not alleles:ctr['nonACGT']+=1;continue
  b=next(iter(alleles));lab='ref' if b==marker['ref'] else 'alt' if b==marker['alt'] else 'other';ctr[lab]+=1
  if any(p>5 for bb,q,m,p in items if bb==b):edge[lab]+=1
 depth=ctr['ref']+ctr['alt'];af=ctr['alt']/depth if depth else None
 return {**marker,'read_ref':read[marker['ref']],'read_alt':read[marker['alt']],'fragment_ref':ctr['ref'],'fragment_alt':ctr['alt'],'fragment_other':ctr['other'],'fragment_discordant':ctr['discordant'],'fragment_nonACGT':ctr['nonACGT'],'fragment_depth_refalt':depth,'fragment_ALT_fraction':af,'minor_allele_fraction':min(af,1-af) if af is not None else None,'alt_forward_reads':strand[(marker['alt'],'fwd')],'alt_reverse_reads':strand[(marker['alt'],'rev')],'ref_forward_reads':strand[(marker['ref'],'fwd')],'ref_reverse_reads':strand[(marker['ref'],'rev')],'fragment_ref_BP5_gt5':edge['ref'],'fragment_alt_BP5_gt5':edge['alt'],'reported_pileup_depth':int(a[3])}
def smoke():
 d=OUT/'synthetic-control';d.mkdir(exist_ok=True);header={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chrTest','LN':1000}],'RG':[{'ID':'test','SM':'test'}]};reads=[]
 for name,alleles,extra,mq,bq in [('ref','AA',0,60,35),('alt','GG',0,60,35),('conflict','AG',0,60,35),('dup','GG',1024,60,35),('supp','GG',2048,60,35),('lowmap','GG',0,20,35),('lowbase','GG',0,60,10)]:
  for j,base in enumerate(alleles):
   r=pysam.AlignedSegment();r.query_name=name;r.flag=(99 if j==0 else 147)|extra;r.reference_id=0;r.reference_start=70;r.mapping_quality=mq;r.cigar=[(0,60)];seq=list('A'*60);seq[30]=base;r.query_sequence=''.join(seq);q=[35]*60;q[30]=bq;r.query_qualities=q;r.next_reference_id=0;r.next_reference_start=70;r.template_length=60 if j==0 else -60;r.set_tag('RG','test');reads.append(r)
 with pysam.AlignmentFile(str(d/'control.bam'),'wb',header=header) as f:
  for r in reads:f.write(r)
 pysam.index(str(d/'control.bam'));(d/'positions.bed').write_text('chrTest\t100\t101\n');args=command(d/'control.bam',d/'positions.bed',d/'control.pileup');pysam.mpileup(*args,catch_stdout=False)
 ln=(d/'control.pileup').read_text();res=parse_line(ln,{'chrom':'chrTest','pos1':101,'ref':'A','alt':'G'})
 zero=parse_line('chrTest\t102\tN\t0\t*\t*\t*\t*\t*\n',{'chrom':'chrTest','pos1':102,'ref':'A','alt':'G'});assert zero['fragment_depth_refalt']==0 and zero['fragment_ALT_fraction'] is None
 assert (res['fragment_ref'],res['fragment_alt'],res['fragment_discordant'])==(1,1,1),res
 assert (res['read_ref'],res['read_alt'])==(3,3),res
 (d/'result.json').write_text(json.dumps({'pysam':pysam.__version__,'samtools':pysam.__samtools_version__,'command':args,'result':res,'status':'PASS','scope':'synthetic filter/mate-collapse/parser execution qualification only, not clinical accuracy'},indent=2)+'\n');print('syntheticcontrolPASS',flush=True)
def run():
 started=time.time();smoke();panel=list(csv.DictReader((OUT/'public-common-SNP-panel.tsv').open(),delimiter='\t'))
 for r in panel:r['pos1']=int(r['pos1'])
 lookup={(r['chrom'],r['pos1']):r for r in panel};summary={}
 for ch in CHROMS:
  dest=OUT/f'common-SNP-alleles.{ch}.tsv'
  if dest.exists():continue
  raw=OUT/f'common-SNP-alleles.{ch}.pileup';args=command(SRC/'DNA_TN26-279853.bam',OUT/f'public-common-SNP-panel.{ch}.bed',raw,ROOT/'work/oct1-analysis/DNA_TN26-279853.bam.bai',ch);t=time.time();pysam.mpileup(*args,catch_stdout=False)
  rows=[]
  with raw.open() as f:
   for ln in f:
    a=ln.split('\t',3);key=(a[0],int(a[1]));rows.append(parse_line(ln,lookup[key]))
  write_tsv(dest,rows);summary[ch]={'selected':sum(r['chrom']==ch for r in panel),'observed_positions':len(rows),'elapsed_seconds':time.time()-t,'command':args,'pileup_bytes':raw.stat().st_size};print('pileup',ch,summary[ch]['observed_positions'],round(time.time()-t,1),flush=True)
  # Preserve compressed raw pileup (query names stay private), remove only this
  # reproducible, newly created intermediate after verifying gzip roundtrip hash.
  h=hashlib.sha256()
  with raw.open('rb') as f,gzip.open(str(raw)+'.gz','wb') as g:
   for z in iter(lambda:f.read(1024*1024),b''):h.update(z);g.write(z)
  hh=hashlib.sha256()
  with gzip.open(str(raw)+'.gz','rb') as f:
   for z in iter(lambda:f.read(1024*1024),b''):hh.update(z)
  assert h.digest()==hh.digest();summary[ch]['raw_sha256']=h.hexdigest();raw.unlink()
  (OUT/'common-SNP-recount.progress.json').write_text(json.dumps(summary,indent=2)+'\n')
 (OUT/'common-SNP-recount.method.json').write_text(json.dumps({'elapsed_seconds':time.time()-started,'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'pysam':pysam.__version__,'samtools':pysam.__samtools_version__,'baseQ_min':25,'MAPQ_min':30,'BAQ':'disabled','overlap_handling':'samtools overlap removal disabled; independently collapse by query name, discard fragments with conflicting ACGT bases','max_depth':100000,'qualifier':'query-name fragments, not UMI unique molecules','read_position':'5-prime distance recorded; >5bp sensitivity recorded; 3-prime-end sensitivity handled by independent direct recount subset','summary':summary},indent=2)+'\n')
if __name__=='__main__':
 if sys.argv[1]=='select':select()
 elif sys.argv[1]=='smoke':smoke()
 elif sys.argv[1]=='run':run()
