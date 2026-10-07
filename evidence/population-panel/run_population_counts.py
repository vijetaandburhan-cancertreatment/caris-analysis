import pathlib,json,csv,hashlib,sys,time,gzip,collections,resource,shutil,subprocess,importlib.util
import pysam
OUT=pathlib.Path(__file__).resolve().parent
WS=pathlib.Path('/Users/burhanazeem/Documents/Codex/2026-09-05/finances-plugin-finances-openai-curated-remote-3')
BASE=OUT.parent
SRC=BASE/'TN26-279853'
CHROMS=['chr2','chr3','chr5','chr8','chr9','chr17']
spec=importlib.util.spec_from_file_location('counter',OUT/'public_counter_source.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def js(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
def rows():
 a=list(csv.DictReader((OUT/'input/public-common-SNP-panel.tsv').open(),delimiter='\t'))
 for r in a:r['pos1']=int(r['pos1'])
 assert len(a)==165782 and len({(r['chrom'],r['pos1']) for r in a})==165782
 return a
def args(bam,bed,raw,index=None,ch=None):
 a=mod.command(bam,bed,raw,index,ch);a[a.index('-d')+1]='0';return a
def qualify():
 control=OUT/'qualification';control.mkdir(exist_ok=True);mod.OUT=control;mod.smoke()
 h={'HD':{'VN':'1.6','SO':'coordinate'},'SQ':[{'SN':'chrTest','LN':1000}],'RG':[{'ID':'test','SM':'test'}]}
 bam=control/'depth-control.bam';bed=control/'depth-control.bed';bed.write_text('chrTest\t100\t101\n')
 with pysam.AlignmentFile(str(bam),'wb',header=h) as f:
  for i in range(9001):
   for j in range(2):
    r=pysam.AlignedSegment();r.query_name=f'control{i}';r.flag=99 if j==0 else 147;r.reference_id=0;r.reference_start=70;r.mapping_quality=60;r.cigar=[(0,60)];r.query_sequence='A'*60;r.query_qualities=[35]*60;r.next_reference_id=0;r.next_reference_start=70;r.template_length=60 if j==0 else -60;r.set_tag('RG','test');f.write(r)
 pysam.index(str(bam));res={}
 for cap in [100,0]:
  raw=control/f'depth-{cap}.pileup';a=args(bam,bed,raw);a[a.index('-d')+1]=str(cap);pysam.mpileup(*a,catch_stdout=False)
  r=mod.parse_line(raw.read_text(),{'chrom':'chrTest','pos1':101,'ref':'A','alt':'G'});res[str(cap)]={'command':a,'result':r}
 assert res['0']['result']['fragment_ref']==9001 and res['0']['result']['reported_pileup_depth']==18002
 assert res['100']['result']['fragment_ref']<9001
 ref=pysam.FastaFile(str(BASE/'genome-fusion-reference/GRCh38.primary_assembly.genome.fa'))
 panel=rows();mismatch=[]
 for r in panel:
  if ref.fetch(r['chrom'],r['pos1']-1,r['pos1']).upper()!=r['ref']:mismatch.append(r)
 assert not mismatch
 bamheaders={}
 for assay in ['DNA','RNA']:
  with pysam.AlignmentFile(str(SRC/f'{assay}_TN26-279853.bam'),'rb',index_filename=str(WS/f'work/oct1-analysis/{assay}_TN26-279853.bam.bai')) as b:
   assert len(b.header.to_dict()['RG'])==1
   assert all(b.get_reference_length(ch)==ref.get_reference_length(ch) for ch in CHROMS)
   bamheaders[assay]={'RG':b.header.to_dict()['RG'],'contig_lengths':{ch:b.get_reference_length(ch) for ch in CHROMS},'header_sha256':hashlib.sha256(str(b.header).encode()).hexdigest(),'bytes':(SRC/f'{assay}_TN26-279853.bam').stat().st_size}
 result={'status':'PASS','samtools':pysam.__samtools_version__,'pysam':pysam.__version__,'public_method_source':'https://www.htslib.org/doc/samtools-mpileup.html','semantics':'Samtools1.24 -d0 sets highest possible depth; 9001 synthetic pairs /18002 records retained. Finite100 cap drops records.','depth_control':res,'reference_checked_loci':len(panel),'reference_mismatches':mismatch,'bam_headers':bamheaders,'completed_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
 js(control/'qualification.json',result)
 d=json.loads((OUT/'design-proposed.json').read_text());d['version']='design-v1-frozen-before-patient-recount';d['frozen_utc']=result['completed_utc'];d['primary_filters']['max_depth']='0; uncapped after synthetic and samtools1.24 manual qualification';d['primary_filters'].pop('reuse_DNA');d['primary_filters']['fresh_DNA']='Recount all public loci with uncapped DNA and RNA. Compare all prior observed DNA counters; prior missing rows remain in universe.';d['callability']['exclude_depth_cap']='Uncapped execution validated; no finite-depth truncation assumption.';d['operational_feasibility_gate'].pop('min_gene_labels');d['operational_feasibility_gate']['min_uniquely_assigned_GENCODE_gene_IDs']=30;d['qualification_sha256']=sha(control/'qualification.json');d['gene_annotation_sha256']=sha(OUT/'public-annotation/public-panel-gene-annotations.tsv');d['public_method_source']=result['public_method_source'];d['run_script_sha256']=sha(pathlib.Path(__file__));d['spacing_rule']='Sort all165782 public sites by numeric chromosome then pos1, greedily retain each if >=50000bp beyond prior retained site. Frozen before patient counts; report callability/dominance in this full-universe subset.'
 js(OUT/'design-frozen.json',d);(OUT/'design-frozen.sha256').write_text(sha(OUT/'design-frozen.json')+'  design-frozen.json\n')
 print('Qualification PASS; design frozen',sha(OUT/'design-frozen.json'),flush=True)
def worker(assay,ch):
 assert (OUT/'design-frozen.json').exists()
 folder=OUT/assay;folder.mkdir(exist_ok=True);panel=[r for r in rows() if r['chrom']==ch];lookup={r['pos1']:r for r in panel}
 bed=folder/f'{ch}.bed';bed.write_text(''.join(f"{ch}\t{r['pos1']-1}\t{r['pos1']}\n" for r in panel))
 raw=folder/f'{ch}.pileup';dest=folder/f'{ch}.tsv';a=args(SRC/f'{assay}_TN26-279853.bam',bed,raw,WS/f'work/oct1-analysis/{assay}_TN26-279853.bam.bai',ch)
 t=time.time();pysam.mpileup(*a,catch_stdout=False);result=[]
 with raw.open() as f:
  for ln in f:
   z=ln.split('\t',3);assert z[0]==ch;r=mod.parse_line(ln,lookup[int(z[1])]);result.append(r)
 mod.write_tsv(dest,result)
 h=hashlib.sha256()
 with raw.open('rb') as f,gzip.open(str(raw)+'.gz','wb',compresslevel=3) as g:
  for z in iter(lambda:f.read(1048576),b''):h.update(z);g.write(z)
 hh=hashlib.sha256()
 with gzip.open(str(raw)+'.gz','rb') as f:
  for z in iter(lambda:f.read(1048576),b''):hh.update(z)
 assert h.digest()==hh.digest();rawbytes=raw.stat().st_size;raw.unlink()
 rec={'assay':assay,'chrom':ch,'command':a,'selected_loci':len(panel),'observed_rows':len(result),'elapsed_seconds':time.time()-t,'maxrss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'raw_sha256':h.hexdigest(),'raw_bytes':rawbytes,'gzip_sha256':sha(str(raw)+'.gz'),'tsv_sha256':sha(dest),'design_sha256':sha(OUT/'design-frozen.json'),'status':'PASS'}
 js(folder/f'{ch}.run.json',rec);print(assay,ch,len(result),round(time.time()-t,1),flush=True)
def run():
 assert sha(OUT/'design-frozen.json')==(OUT/'design-frozen.sha256').read_text().split()[0]
 progress=[]
 for assay in ['DNA','RNA']:
  for ch in CHROMS:
   complete=OUT/assay/f'{ch}.run.json'
   if complete.exists():progress.append(json.loads(complete.read_text()));continue
   before=shutil.disk_usage(OUT).free;assert before>3*1024**3,'Disk guard before count'
   log=OUT/f'{assay}-{ch}.log';t=time.time()
   with log.open('w') as f:
    p=subprocess.Popen([sys.executable,__file__,'worker',assay,ch],stdout=f,stderr=subprocess.STDOUT)
    while p.poll() is None:
     time.sleep(1)
     rss=subprocess.run(['ps','-o','rss=','-p',str(p.pid)],capture_output=True,text=True).stdout.strip();rss=int(rss or 0)*1024
     free=shutil.disk_usage(OUT).free
     artifacts=sum(z.stat().st_size for z in OUT.rglob('*') if z.is_file())
     if rss>4*1024**3 or free<3*1024**3 or time.time()-t>1800 or artifacts>1024**3:
      p.terminate();p.wait();raise RuntimeError(f'Guard triggered RSS={rss}, free={free}, artifacts={artifacts}, seconds={time.time()-t}')
   assert p.returncode==0,log.read_text();progress.append(json.loads(complete.read_text()));js(OUT/'count-progress.json',progress);print(log.read_text().strip(),flush=True)
 js(OUT/'count-complete.json',{'status':'PASS','runs':len(progress),'seconds':sum(x['elapsed_seconds'] for x in progress),'design_sha256':sha(OUT/'design-frozen.json'),'runs_detail':progress})
if __name__=='__main__':
 if sys.argv[1]=='qualify':qualify()
 elif sys.argv[1]=='run':run()
 elif sys.argv[1]=='worker':worker(sys.argv[2],sys.argv[3])
