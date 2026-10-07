#!/usr/bin/env python3
"""Build full public G37 STAR sparse index with bounded memory and disk guard."""
import pathlib,subprocess,time,json,os,shutil,signal,hashlib,datetime,sys
WORK=pathlib.Path(__file__).resolve().parent
CACHE=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis')
REF=CACHE/'genome-fusion-reference'; IDX=REF/'STAR_G37_primary_D8_SA12'; IDX.mkdir(exist_ok=True)
STAR=CACHE/'genome-fusion-tools/STAR_2.7.11b/MacOSX_x86_64/STAR'
if (IDX/'complete.json').exists():print('Index already marked complete');sys.exit(0)
assert not (IDX/'Genome').exists(),'Refuse ambiguous partial/previous index overwrite; inspect first'
args=[str(STAR),'--runMode','genomeGenerate','--runThreadN','4','--genomeDir',str(IDX),'--genomeFastaFiles',str(REF/'GRCh38.primary_assembly.genome.fa'),'--sjdbGTFfile',str(REF/'gencode.v37.primary_assembly.annotation.gtf'),'--sjdbOverhang','160','--genomeSAsparseD','8','--genomeSAindexNbases','12','--limitGenomeGenerateRAM','22000000000','--outFileNamePrefix',str(IDX/'build.')]
assert shutil.disk_usage(REF).free>11*1024**3, 'Need >11GiB free to start full index safely'
state={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stage':'indexing','command':args,'index':str(IDX),'disk_floor_bytes':3*1024**3,'RSS_guard_bytes':27*1024**3,'STAR_sha256':hashlib.sha256(STAR.read_bytes()).hexdigest(),'STAR_version':subprocess.check_output([str(STAR),'--version'],text=True).strip(),'full_primary_assembly':True,'raw_patient_files_read':False,'algorithm_source_modified':False}
log=(IDX/'driver.log').open('w');p=subprocess.Popen(args,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);state['pid']=p.pid
caff=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())]);begin=time.time();peak=0
try:
 while p.poll() is None:
  free=shutil.disk_usage(REF).free
  try:rss=int(subprocess.check_output(['ps','-o','rss=','-p',str(p.pid)],text=True).strip())*1024
  except Exception:rss=0
  peak=max(peak,rss);state.update(elapsed_seconds=round(time.time()-begin,1),free_disk_bytes=free,peak_RSS_bytes=peak,updated_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
  if free<state['disk_floor_bytes'] or rss>state['RSS_guard_bytes']:
   state['stage']='guard_stop';state['stop_reason']='disk reserve below3GiB' if free<state['disk_floor_bytes'] else 'STAR RSS above27GiB';os.killpg(p.pid,signal.SIGTERM);p.wait(timeout=20);break
  (WORK/'index-status.json').write_text(json.dumps(state,indent=2));time.sleep(20)
 rc=p.wait();state.update(exit_code=rc,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_seconds=round(time.time()-begin,1),peak_RSS_bytes=peak)
 if rc==0:
  assert all((IDX/f).stat().st_size>0 for f in ['Genome','SA','SAindex','genomeParameters.txt','sjdbInfo.txt'])
  state['stage']='complete';state['index_file_sizes']={f.name:f.stat().st_size for f in IDX.iterdir() if f.is_file()};(IDX/'complete.json').write_text(json.dumps(state,indent=2))
 elif state['stage']!='guard_stop':state['stage']='failed'
except BaseException as e:
 state['exception']=repr(e);state['stage']='failed'
 if p.poll() is None:
  os.killpg(p.pid,signal.SIGTERM)
  try:p.wait(timeout=10)
  except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL)
 raise
finally:
 (WORK/'index-status.json').write_text(json.dumps(state,indent=2));log.close();caff.terminate()
print(json.dumps(state,indent=2))
assert state['stage']=='complete',state
