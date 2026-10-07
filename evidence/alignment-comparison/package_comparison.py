from pathlib import Path
import json,hashlib,zipfile,datetime,shutil
R=Path(__file__).resolve().parent;D=R/'delivery';D.mkdir(exist_ok=True)
assert json.loads((R/'patient-full-pass/run.json').read_text())['status']=='COMPLETE'
assert json.loads((D/'final-content-approval.json').read_text())['approved'] is True
# Immutable scientific evidence package. AWS/account administration stays separate.
files=[]
for sub in ['capture-qualification','independent-feasibility','comparison','name-comparison','independent-comparison','original-complete-names','design-history']:
 files += [p for p in (R/sub).rglob('*') if p.is_file() and '__pycache__' not in p.parts and not p.name.endswith('.pyc')]
files += [p for p in R.glob('*') if p.is_file() and p.suffix in ['.py','.json','.tsv','.bed','.txt','.sha256','.log']]
for name in ['run.json','capture-summary.json','cohort-and-exception-records.bam','STAR.Log.final.out','STAR.Log.out','STAR.Log.progress.out','capture.log','STAR.stderr.log']:
 files.append(R/'patient-full-pass'/name)
files += [p for p in D.glob('*') if p.is_file() and p.suffix in ['.pdf','.txt','.json','.py'] and '-preview' not in p.name and p.name not in ['package-manifest.json','package-verification.json','excluded-resident-files.json']]
files=sorted(set(files));assert all(p.exists() for p in files)
def digest(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
all_run=[p for p in (R/'patient-full-pass').rglob('*') if p.is_file()]
excluded=[{'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':digest(p),'reason':'Not required for allele/placement comparison; preserved resident, not deleted.'} for p in all_run if p not in files]
(D/'excluded-resident-files.json').write_text(json.dumps(excluded,indent=2)+'\n');files.append(D/'excluded-resident-files.json')
manifest={'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Final fixed-cohort full-pair alignment comparison. Original Caris data and reference remain resident; previous milestone packages remain immutable.','files':[{'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':digest(p)} for p in files]}
mp=D/'package-manifest.json';mp.write_text(json.dumps(manifest,indent=2)+'\n');files.append(mp)
assert shutil.disk_usage(D).free>3.5*1024**3
z=D/'Caris-final-alignment-comparison-evidence-2026-10-05.zip';assert not z.exists()
with zipfile.ZipFile(z,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as a:
 for p in files:a.write(p,str(p.relative_to(R)),compress_type=zipfile.ZIP_STORED if p.suffix in ['.bam','.gz','.pdf'] else zipfile.ZIP_DEFLATED)
with zipfile.ZipFile(z) as a:
 assert a.testzip() is None
 for x in manifest['files']:
  data=a.read(x['path']);assert len(data)==x['bytes'];assert hashlib.sha256(data).hexdigest()==x['sha256']
 assert a.read('delivery/package-manifest.json')==mp.read_bytes()
receipt={'status':'PASS','zip':str(z),'bytes':z.stat().st_size,'sha256':digest(z),'archive_members':len(files),'all_member_SHA256_verified':True,'manifest_sha256':digest(mp),'excluded_resident_files':len(excluded),'remaining_free_bytes':shutil.disk_usage(D).free};(D/'package-verification.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
