#!/usr/bin/env python3
"""Default Arriba-recommended STAR alignment, fully streamed with resource guards."""
import pathlib,subprocess,time,json,os,shutil,signal,hashlib,datetime,argparse,sys
CACHE=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis');TOOLS=CACHE/'genome-fusion-tools';STAR=TOOLS/'STAR_2.7.11b/MacOSX_x86_64/STAR';ARRIBA=TOOLS/'arriba_v2.5.1/arriba';DB=TOOLS/'arriba_v2.5.1/database'
COMMON=['--genomeLoad','NoSharedMemory','--readFilesCommand','/usr/bin/gzip','-cd','--outStd','BAM_Unsorted','--outSAMtype','BAM','Unsorted','--outSAMunmapped','Within','--outBAMcompression','0','--outFilterMultimapNmax','50','--peOverlapNbasesMin','10','--alignSplicedMateMapLminOverLmate','0.5','--alignSJstitchMismatchNmax','5','-1','5','5','--chimSegmentMin','10','--chimOutType','WithinBAM','HardClip','Junctions','--chimJunctionOverhangMin','10','--chimScoreDropMax','30','--chimScoreJunctionNonGTAG','0','--chimScoreSeparation','1','--chimSegmentReadGapMax','3','--chimMultimapNmax','50','--chimOutJunctionFormat','1']

def terminate(procs):
 for p in procs:
  if p.poll() is None:
   try:os.killpg(p.pid,signal.SIGTERM)
   except ProcessLookupError:pass
 for p in procs:
  try:p.wait(timeout=10)
  except subprocess.TimeoutExpired:
   try:os.killpg(p.pid,signal.SIGKILL)
   except ProcessLookupError:pass

def run(index,reference,gtf,r1,r2,out,threads=6,save_bam=False,expected_pairs=None,capture_bed=None):
 out=pathlib.Path(out);out.mkdir(parents=True,exist_ok=True)
 assert not (out/'run.json').exists(),'Refuse overwriting previous run; use a new directory'
 args=[str(STAR),'--runThreadN',str(threads),'--genomeDir',str(index),'--readFilesIn',str(r1),str(r2),'--outFileNamePrefix',str(out/'STAR.')]+COMMON
 acmd=[str(ARRIBA),'-x','/dev/stdin','-o',str(out/'fusions.tsv'),'-O',str(out/'fusions.discarded.tsv'),'-a',str(reference),'-g',str(gtf),'-b',str(DB/'blacklist_hg38_GRCh38_v2.5.1.tsv.gz'),'-k',str(DB/'known_fusions_hg38_GRCh38_v2.5.1.tsv.gz'),'-t',str(DB/'known_fusions_hg38_GRCh38_v2.5.1.tsv.gz'),'-p',str(DB/'protein_domains_hg38_GRCh38_v2.5.1.gff3')]
 meta={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'running','STAR_command':args,'Arriba_command':acmd,'save_full_BAM':save_bam,'caller_filters':'All Arriba defaults; supplied blacklist, known-fusion recovery/tags and domains. No -f disabling filters, no threshold changes.','STAR_extra_vs_recommended':'Junctions plus chimOutJunctionFormat1 for audit output; gzip -cd instead of zcat for macOS; remaining recommended mapping parameters retained.','stage':'STAR streaming all alignments to Arriba','expected_input_pairs':expected_pairs}
 err=(out/'STAR.stderr.log').open('w');alog=(out/'Arriba.log').open('w');procs=[];t0=time.time();peak=0;caff=None;fifo=None;capture_log=None
 try:
  s=subprocess.Popen(args,stdout=subprocess.PIPE,stderr=err,start_new_session=True);procs.append(s);stream=s.stdout
  if capture_bed:
   fifo=out/'target-capture.fifo';os.mkfifo(fifo);capture_log=(out/'target-capture.log').open('w')
   capture_cmd=[sys.executable,str(pathlib.Path(__file__).resolve().parent/'capture-target-bam.py'),str(capture_bed),str(fifo),str(out/'target-genes.bam')]
   reader=subprocess.Popen(capture_cmd,stdout=capture_log,stderr=subprocess.STDOUT,start_new_session=True);procs.append(reader)
   capture_tee=subprocess.Popen(['/usr/bin/tee',str(fifo)],stdin=stream,stdout=subprocess.PIPE,start_new_session=True);stream.close();stream=capture_tee.stdout;procs.append(capture_tee)
   meta['selective_BAM_capture']={'command':capture_cmd,'bed':str(capture_bed),'bed_sha256':hashlib.sha256(pathlib.Path(capture_bed).read_bytes()).hexdigest(),'scope':'All alignment records overlapping supplied genomic intervals; mates outside intervals not guaranteed retained.'}
  if save_bam:
   tee=subprocess.Popen(['/usr/bin/tee',str(out/'Aligned.control.bam')],stdin=stream,stdout=subprocess.PIPE,start_new_session=True);stream.close();stream=tee.stdout;procs.append(tee)
  a=subprocess.Popen(acmd,stdin=stream,stdout=alog,stderr=subprocess.STDOUT,start_new_session=True);stream.close();procs.append(a);caff=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())]);meta['processes']=[{'pid':p.pid} for p in procs]
  while any(p.poll() is None for p in procs):
   free=shutil.disk_usage(out).free
   rss=0
   for p in procs:
    if p.poll() is None:
     try:rss+=int(subprocess.check_output(['ps','-o','rss=','-p',str(p.pid)],text=True).strip())*1024
     except Exception:pass
   peak=max(peak,rss);meta.update(elapsed_seconds=round(time.time()-t0,2),current_RSS_bytes=rss,peak_RSS_bytes=peak,free_disk_bytes=free,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());(out/'status.json').write_text(json.dumps(meta,indent=2))
   if free<3*1024**3 or rss>27*1024**3:raise RuntimeError('Resource guard: disk below3GiB or total pipeline RSS above27GiB')
   if any(p.poll() is not None and p.returncode!=0 for p in procs):raise RuntimeError('Pipeline subprocess failed: '+str([p.poll() for p in procs]))
   time.sleep(10)
  meta['exit_codes']=[p.wait() for p in procs];assert all(c==0 for c in meta['exit_codes'])
  assert (out/'fusions.tsv').exists() and (out/'fusions.discarded.tsv').exists()
  final={}
  for line in (out/'STAR.Log.final.out').read_text().splitlines():
   if '|' in line:k,v=line.split('|',1);final[k.strip()]=v.strip()
  meta['STAR_final']=final
  if expected_pairs is not None:assert int(final['Number of input reads'])==expected_pairs,(final['Number of input reads'],expected_pairs)
  assert (out/'STAR.Chimeric.out.junction').exists();assert (out/'STAR.SJ.out.tab').exists()
  meta['status']='complete'
 except BaseException as e:
  meta['status']='failed';meta['exception']=repr(e);terminate(procs);raise
 finally:
  meta['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();meta['elapsed_seconds']=round(time.time()-t0,2);meta['peak_RSS_bytes']=peak;meta['files']={f.name:f.stat().st_size for f in out.iterdir() if f.is_file()};(out/'run.json').write_text(json.dumps(meta,indent=2));(out/'status.json').write_text(json.dumps(meta,indent=2));err.close();alog.close()
  if caff:caff.terminate()
  if capture_log:capture_log.close()
  if fifo and fifo.exists():fifo.unlink()
 return meta

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--index',required=True);ap.add_argument('--reference',required=True);ap.add_argument('--gtf',required=True);ap.add_argument('--r1',required=True);ap.add_argument('--r2',required=True);ap.add_argument('--out',required=True);ap.add_argument('--threads',type=int,default=6);ap.add_argument('--save-bam',action='store_true');ap.add_argument('--expected-pairs',type=int);ap.add_argument('--capture-bed');v=ap.parse_args();print(json.dumps(run(v.index,v.reference,v.gtf,v.r1,v.r2,v.out,v.threads,v.save_bam,v.expected_pairs,v.capture_bed),indent=2))
