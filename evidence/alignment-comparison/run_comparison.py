from pathlib import Path
import subprocess,sys,json,time,os,signal,shutil,hashlib,datetime,csv,pysam
R=Path(__file__).resolve().parent;B=R.parent;OUT=R/'patient-full-pass'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write(p,o):p.write_text(json.dumps(o,indent=2)+'\n')
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

def prepare():
 old=json.loads((B/'oct4-followup/genome-fusion/patient-D8/run.json').read_text());cmd=old['STAR_command'][:];cmd[cmd.index('--outFileNamePrefix')+1]=str(OUT/'STAR.')
 assert cmd[cmd.index('--runThreadN')+1]=='6'
 free=shutil.disk_usage(R).free;assert free>=4*1024**3,f'Launch requires4GiB free; now{free}'
 q=json.loads((R/'capture-qualification/qualification.json').read_text());assert q['status']=='PASS'
 design={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'FROZEN_PRELAUNCH','question':'Do baseline RNA allele discrepancies persist under existing full-reference STAR2.7.11b/D8 pipeline? Comparison of methods, not correction or attribution to sparsity.','cohort':'All4015 jointly callable original loci, fixed original DNA categories/denominators; retainall13 original exceptions and3 separateMAPQ60flags evenwhen newcoveragezero.','inputs':json.loads((R/'frozen-selection.json').read_text()),'STAR_command':cmd,'expected_pairs':23209264,'capture':'Compressed BAM union of fixed4015 SNP M/=/X overlaps and all emitted records for frozen exceptionQNAME set across everycontig/flag, including mates/secondary/supplementary/unmapped. No capture qualityfilters. Doesnotretain every possiblemapping.','postcomparison':'OriginalDNAfixed; newRNA matchedMAPQ30/BQ25/properpairs/ff0xF0C; collapseQNAME, reject conflictingACGT. Report rawcounts,coverage losses/fixedcohort fractions, zero/missing sites, novel flags; no replacingdenominator withsurvivors. Directstream counting canavoid sorting. Independent implementation auditsmaterialcounts.','complete_mate_context':'Freezeall names whose originalRNA Mblock covers any16exception/mappingflag sites withoutqualityfilter. Preserveallnewemittedrecords. Retain originalallrecords for these exactnames through separateboundedstreamscan ifneeded.','resource_guards':{'preflight_free_bytes_min':4*1024**3,'disk_abort_bytes':int(3.5*1024**3),'absolute_requested_floor_bytes':3*1024**3,'combined_RSS_abort_bytes':12*1024**3,'capture_abort_bytes':256*1024**2,'own_run_bytes_abort':512*1024**2,'poll_seconds':2,'wall_seconds_abort':3600,'deadline_utc':'2026-10-05T10:15:00Z'},'resource_limit_note':'Disk floor monitored with0.5GiB margin; independent hostwrites/swapallocations canoccur betweenpolls. Stopandpreservepartial/negative result uponlimits; no deletionoforiginal/evidence orcloudworkaround.','source_differences':'OriginalSTAR2.7.8a CTAT/twopassBasic/flushRight/overlap12 etc vs local2.7.11b primaryassembly/D8/onepass/localparameters; notsinglefactorcomparison. Capturequalityfiltersmatch, upstreamflagging/mapping maydiffer.','prohibited':'NoArriba, fullBAM, new clinicalmutation/identitycertainty/contaminationfraction/target/sensitivityclaim, externalmessages, AWSupgrade/IAM orreset.','capture_script_sha256':sha(R/'select_capture.py'),'run_script_sha256':sha(Path(__file__)),'capture_qualification_sha256':sha(R/'capture-qualification/qualification.json'),'initial_free_bytes':free}
 write(R/'design-frozen.json',design);(R/'design-frozen.sha256').write_text(sha(R/'design-frozen.json')+'  design-frozen.json\n');print('Design frozen',sha(R/'design-frozen.json'))
def run():
 d=json.loads((R/'design-frozen.json').read_text());assert sha(R/'design-frozen.json')==(R/'design-frozen.sha256').read_text().split()[0];assert sha(R/'select_capture.py')==d['capture_script_sha256'];assert sha(Path(__file__))==d['run_script_sha256'];audit=json.loads((R/'independent-feasibility/launch-audit.json').read_text());assert audit['status']=='APPROVED' and audit['design_sha256']==sha(R/'design-frozen.json')
 for key,filename in [('cohort','baseline-fixed-4015.tsv'),('exception_sites','frozen-exception-sites.tsv'),('exception_qnames','frozen-exception-qnames.txt')]:assert sha(R/filename)==d['inputs']['input_sha256'][key]
 cohort=list(csv.DictReader((R/'baseline-fixed-4015.tsv').open(),delimiter='\t'));expected=[(z['chrom'],int(z['pos1'])-1,int(z['pos1'])) for z in cohort];actual=[(ch,int(a),int(b)) for ch,a,b in (ln.split('\t') for ln in (R/'baseline-sites.bed').read_text().splitlines())];assert actual==expected and len(expected)==4015
 assert not OUT.exists();OUT.mkdir();g=d['resource_guards'];assert shutil.disk_usage(R).free>=g['preflight_free_bytes_min'];procs=[];t=time.time();peak=0;minfree=shutil.disk_usage(R).free
 status={'status':'RUNNING','design_sha256':sha(R/'design-frozen.json'),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};caf=None
 with (OUT/'STAR.stderr.log').open('w') as err,(OUT/'capture.log').open('w') as log:
  try:
   s=subprocess.Popen(d['STAR_command'],stdout=subprocess.PIPE,stderr=err,start_new_session=True);procs.append(s)
   capcmd=[sys.executable,str(R/'select_capture.py'),'-',str(R/'baseline-sites.bed'),str(R/'frozen-exception-qnames.txt'),str(OUT/'cohort-and-exception-records.bam'),str(OUT/'capture-summary.json')]
   c=subprocess.Popen(capcmd,stdin=s.stdout,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);s.stdout.close();procs.append(c);status['capture_command']=capcmd;status['pids']=[p.pid for p in procs];caf=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())])
   while any(p.poll() is None for p in procs):
    free=shutil.disk_usage(R).free;minfree=min(minfree,free);rss=0
    for p in procs:
     if p.poll() is None:
      v=subprocess.run(['ps','-o','rss=','-p',str(p.pid)],capture_output=True,text=True).stdout.strip();rss+=int(v or 0)*1024
    peak=max(peak,rss);total=sum(f.stat().st_size for f in OUT.rglob('*') if f.is_file());cap=OUT/'cohort-and-exception-records.bam';cb=cap.stat().st_size if cap.exists() else 0
    status.update(elapsed_seconds=time.time()-t,current_RSS_bytes=rss,peak_RSS_bytes=peak,free_disk_bytes=free,min_free_disk_bytes=minfree,own_run_bytes=total,capture_bytes=cb,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());write(OUT/'status.json',status)
    if free<g['disk_abort_bytes'] or rss>g['combined_RSS_abort_bytes'] or cb>g['capture_abort_bytes'] or total>g['own_run_bytes_abort'] or time.time()-t>g['wall_seconds_abort'] or datetime.datetime.now(datetime.timezone.utc)>=datetime.datetime.fromisoformat(g['deadline_utc'].replace('Z','+00:00')):raise RuntimeError('Resource guard: '+json.dumps(status))
    if any(p.poll() is not None and p.returncode!=0 for p in procs):raise RuntimeError('Pipeline nonzeroexit '+str([p.poll() for p in procs]))
    time.sleep(g['poll_seconds'])
   assert all(p.wait()==0 for p in procs)
   free=shutil.disk_usage(R).free;minfree=min(minfree,free);total=sum(f.stat().st_size for f in OUT.rglob('*') if f.is_file());cb=(OUT/'cohort-and-exception-records.bam').stat().st_size
   assert free>=g['disk_abort_bytes'] and cb<=g['capture_abort_bytes'] and total<=g['own_run_bytes_abort'],'Final-close resource guard'
   pysam.quickcheck(str(OUT/'cohort-and-exception-records.bam'));status['final_closed_BAM_quickcheck']='PASS';status['final_close_free_disk_bytes']=free;status['final_closed_capture_bytes']=cb;sf={}
   for ln in (OUT/'STAR.Log.final.out').read_text().splitlines():
    if '|' in ln:k,v=ln.split('|',1);sf[k.strip()]=v.strip()
   assert int(sf['Number of input reads'])==d['expected_pairs'];capture=json.loads((OUT/'capture-summary.json').read_text());assert capture['status']=='CAPTURE_STREAM_EOF';status.update(status='COMPLETE',STAR_final=sf,capture_summary=capture,exit_codes=[p.returncode for p in procs])
  except BaseException as e:
   status.update(status='STOPPED_OR_FAILED',exception=repr(e));terminate(procs);raise
  finally:
   status.update(finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_seconds=time.time()-t,peak_RSS_bytes=peak,min_free_disk_bytes=minfree);write(OUT/'run.json',status);write(OUT/'status.json',status)
   if caf:caf.terminate()
 print('Fullpass',status['status'],round(status['elapsed_seconds'],1),flush=True)
if __name__=='__main__':
 if sys.argv[1]=='prepare':prepare()
 elif sys.argv[1]=='run':run()
