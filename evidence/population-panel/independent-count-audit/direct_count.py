"""Independent indexed CIGAR-block allele audit. No pileup/shared parser imported."""
from pathlib import Path
import csv,json,bisect,hashlib,collections,time,resource,sys,datetime
import pysam
OUT=Path(__file__).resolve().parent;ROOT=OUT.parent
WS=Path('/Users/burhanazeem/Documents/Codex/2026-09-05/finances-plugin-finances-openai-curated-remote-3')
DATA=ROOT.parent/'TN26-279853'
EXPECTED_DESIGN='90a153f630375103363792037618d93a864b01bad029575f36d6a3a382ee0bfa'
MODES=['baseline','clean_no_softclip_edge5','MAPQ60','NH_unique_if_tagged']
MASK={b:1<<i for i,b in enumerate('ACGT')}
CHR=['chr2','chr3','chr5','chr8','chr9','chr17']
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def records(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def js(p,o):p.write_text(json.dumps(o,indent=2)+'\n')
def windows(sites,gap=500,maxspan=50000):
 out=[]
 for r in sites:
  if out and r['pos1']-out[-1][-1]['pos1']<=gap and r['pos1']-out[-1][0]['pos1']<maxspan:out[-1].append(r)
  else:out.append([r])
 return out

def count_window(bam,chrom,sites):
 pos=[r['pos1']-1 for r in sites]
 groups=[[{} for p in pos] for mode in MODES]
 reads=[collections.Counter() for p in pos]
 inspected=eligible=0
 for rd in bam.fetch(chrom,pos[0],pos[-1]+1):
  inspected+=1
  if rd.flag&0xF0C or not rd.flag&2 or rd.mapping_quality<30:continue
  seq=rd.query_sequence;qual=rd.query_qualities
  if seq is None or qual is None:continue
  eligible+=1
  cigar=rd.cigartuples or []
  no_soft=not any(op==4 for op,n in cigar)
  unique=not rd.has_tag('NH') or rd.get_tag('NH')==1
  rp=rd.reference_start;qp=0
  for op,n in cigar:
   if op in (0,7,8):
    lo=bisect.bisect_left(pos,rp);hi=bisect.bisect_left(pos,rp+n)
    for idx in range(lo,hi):
     q=qp+pos[idx]-rp
     if qual[q]<25:continue
     base=seq[q].upper()
     if base not in MASK:continue
     bit=MASK[base];name=rd.query_name
     reads[idx][base]+=1
     passes=[True,no_soft and q-rd.query_alignment_start>=5 and rd.query_alignment_end-1-q>=5,rd.mapping_quality>=60,unique]
     for j,yes in enumerate(passes):
      if yes:
       g=groups[j][idx];g[name]=g.get(name,0)|bit
    rp+=n;qp+=n
   elif op in (1,4):qp+=n
   elif op in (2,3):rp+=n
   elif op in (5,6):pass
   else:raise RuntimeError(f'Unhandled CIGAR op{op}')
 results=[]
 for i,r in enumerate(sites):
  z={'chrom':chrom,'pos1':r['pos1'],'id':r.get('id','.'),'ref':r['ref'],'alt':r['alt']}
  for mode,gs in zip(MODES,groups):
   ct=collections.Counter(gs[i].values());ctbase={b:ct[mask] for b,mask in MASK.items()}
   z.update({f'{mode}_{b}':ctbase[b] for b in 'ACGT'})
   z[f'{mode}_fragment_ref']=ctbase[r['ref']];z[f'{mode}_fragment_alt']=ctbase[r['alt']]
   z[f'{mode}_fragment_other']=sum(ctbase[b] for b in 'ACGT' if b not in (r['ref'],r['alt']))
   z[f'{mode}_fragment_discordant']=sum(v for m,v in ct.items() if m not in MASK.values())
   z[f'{mode}_ACGT_depth']=sum(ctbase.values())
  z['read_ref']=reads[i][r['ref']];z['read_alt']=reads[i][r['alt']]
  results.append(z)
 return results,{'fetched_records':inspected,'qualified_read_records':eligible}

def run(assay):
 assert (ROOT/'count-complete.json').exists() and (ROOT/'joint-callable.tsv').exists(),'Wait for primary completion and merged table.'
 assert sha(ROOT/'design-frozen.json')==EXPECTED_DESIGN
 sites=records(ROOT/'joint-callable.tsv')
 for r in sites:r['pos1']=int(r['pos1'])
 assert len({(r['chrom'],r['pos1']) for r in sites})==len(sites)
 sites.sort(key=lambda r:(CHR.index(r['chrom']),r['pos1']))
 provenance={'assay':assay,'scope':'All primary joint-callable positions; direct independent CIGAR-block/count implementation.','design_sha256':EXPECTED_DESIGN,'joint_table_sha256':sha(ROOT/'joint-callable.tsv'),'script_sha256':sha(Path(__file__)),'pysam':pysam.__version__,'samtools':pysam.__samtools_version__,'filters':'MAPQ>=30,BQ>=25,proper flag2,exclude0xF0C,qualified ACGT observations collapsed per QNAME; conflicting ACGT names excluded. Query sequence in BAM alignment orientation; no extra reverse complement.','sensitivity':'Clean requires no S anywhere in read CIGAR and >=5 query bases between locus and each aligned query end. MAPQ60 and NH==1-if-present are separate sensitivity masks. None certifies uniqueness.','excluded_from_exact_comparison':'nonACGT symbolic deletion/refskip accounting is not reconstructed; exact comparison covers ACGT alleles, conflicts and REF/ALT read counts only.'}
 t=time.time();progress=[];total=collections.Counter();dest=OUT/f'{assay}-direct-counts.tsv'
 with pysam.AlignmentFile(str(DATA/f'{assay}_TN26-279853.bam'),'rb',index_filename=str(WS/f'work/oct1-analysis/{assay}_TN26-279853.bam.bai')) as bam,dest.open('w') as f:
  writer=None
  for chrom in CHR:
   cc=[r for r in sites if r['chrom']==chrom];ww=windows(cc);tc=time.time()
   for wi,w in enumerate(ww):
    zz,counts=count_window(bam,chrom,w);total.update(counts);total['sites']+=len(zz)
    if writer is None:writer=csv.DictWriter(f,fieldnames=list(zz[0]),delimiter='\t',lineterminator='\n');writer.writeheader()
    writer.writerows(zz)
    if wi%50==0:
     f.flush();js(OUT/f'{assay}-progress.json',{'chrom':chrom,'window':wi,'windows':len(ww),'sites':total['sites'],'seconds':time.time()-t,'RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
     if time.time()-t>1800:raise RuntimeError('30min audit guard; partial preserved')
     if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>4*1024**3:raise RuntimeError('4GiB RSS guard')
   progress.append({'chrom':chrom,'sites':len(cc),'windows':len(ww),'seconds':time.time()-tc});print(assay,chrom,len(cc),round(time.time()-tc,2),flush=True)
 assert total['sites']==len(sites)
 # Compare separately loaded primary tables rather than importing their parser.
 primary={}
 for chrom in CHR:
  for r in records(ROOT/assay/f'{chrom}.tsv'):primary[(r['chrom'],int(r['pos1']))]=r
 dif=[]
 for r in records(dest):
  p=primary[(r['chrom'],int(r['pos1']))];delta={}
  for field in ('fragment_ref','fragment_alt','fragment_other','fragment_discordant'):
   if int(r['baseline_'+field])!=int(p[field]):delta[field]={'direct':int(r['baseline_'+field]),'primary':int(p[field])}
  for field in ('read_ref','read_alt'):
   if int(r[field])!=int(p[field]):delta[field]={'direct':int(r[field]),'primary':int(p[field])}
  if delta:dif.append({'chrom':r['chrom'],'pos1':int(r['pos1']),'difference':delta})
 result={**provenance,'status':'PASS' if not dif else 'REVIEW_REQUIRED','sites':len(sites),'direct_output_sha256':sha(dest),'exact_compared_fields':['fragment_ref','fragment_alt','fragment_other','fragment_discordant','read_ref','read_alt'],'differences':dif,'progress':progress,'totals':dict(total),'seconds':time.time()-t,'peak_RSS_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 js(OUT/f'{assay}-audit.json',result);print(assay,result['status'],len(dif),'differences',flush=True)
if __name__=='__main__':run(sys.argv[1])
