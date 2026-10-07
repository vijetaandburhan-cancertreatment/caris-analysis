from pathlib import Path
import json,hashlib,zipfile,datetime,shutil

R=Path.cwd();W=Path(__file__).parent;O=R/'outputs/caris-followup';L=Path.home()/'.local/share/codex/caris-analysis'
old=L/'oct5-final-handoff/Caris-computational-followups-2026-10-05.zip'
assert hashlib.sha256(old.read_bytes()).hexdigest()=='143936a401fa47395a30b6574c854212c1f585a5a7f44a8db834f004b72e6a18'
entries={};excluded=[]
def add(p,arc):
 assert p.is_file(),p
 entries[arc]=p

for p in [R/'output/pdf/Caris-DNA-RNA-research-handoff-TN26-279853-2026-10-05-v2.pdf', O/'Caris-followup-findings-2026-10-05-v2.pdf',O/'Caris-followup-findings-2026-10-05-v2.txt',O/'Caris-morning-update-2026-10-05-v2.txt',O/'Caris-specialist-questions-2026-10-05.txt']:
 add(p,p.name)
for branch,arc in [('work/oct5-CCND1-HEG1-audit-v1','CCND1-HEG1-additional'),('work/oct5-fusion-robustness-v1','fusion-robustness')]:
 root=L/Path(branch).name
 assert root.is_dir(),root
 for p in sorted(root.rglob('*')):
  if not p.is_file():continue
  rel=p.relative_to(root)
  if any(x in ('__pycache__','parser-history','audit-parser-history','analysis-cache','analysis-cache-initial-line-parser') for x in rel.parts) or 'partial' in p.name or p.stat().st_size>20_000_000:
   excluded.append({'path':str(rel),'reason':'Intermediate/cache or over20MB; retained separately','bytes':p.stat().st_size});continue
  if p.suffix.lower() in ('.py','.json','.tsv','.csv','.txt','.fa','.fasta','.bam','.gz','.out','.log','.junction','.tab','.sha256'):
   add(p,arc+'/'+str(rel))
for p in W.rglob('*'):
 if p.is_file() and p.suffix in ('.py','.json') and '__pycache__' not in p.parts:add(p,'report-build-v2/'+str(p.relative_to(W)))

readme='''CARIS COMPUTATIONAL FOLLOW-UP EVIDENCE — 5 OCTOBER2026, V2

This version retains the original four-workstream evidence and adds the prespecified45-run public fusion benchmark and bounded CCND1::HEG1 RNA/DNA audits. Start with the revised one-page handoff and four-page report. Earlier delivered artifacts remain unchanged outside this new bundle.

The public benchmark uses tiny overlapping/modified subsets of a known BCR::ABL1 fixture with a small synthetic normal-parent background. It cannot estimate patient-level clinical sensitivity, specificity, detection limit, or abundance. The unchanged frozen design and prospective/technical amendments are included. Stock Arriba no-chimeric-input exit1 is retained as no call table, never fabricated as an empty result.

CCND1::HEG1 sequence support is distinguished from its biological origin. RNA microhomology and highly overlapping mates limit breakpoint and molecule claims; locus-limited DNA checks cannot exclude intronic rearrangements. See the DNA addendum for exact-marker and split-alignment findings. Independent root RNA scripts named affine_dp.py and affine_terminal_dp.py are exploratory end-convention comparisons; affine_window_dp.py reproduces the frozen physical-window convention used for patient scoring. The audit findings describe these distinctions.

This is a compact research handoff, not a clinical assay or a full raw-data archive. Nine original Caris files, references, software/models and larger intermediates remain locally preserved; the original files also have previously verified Drive copies. Every included member except the manifest itself is SHA256-listed in MANIFEST.json. Excluded intermediates are listed in PACKAGE-EXCLUSIONS.json. File paths in provenance refer to the analysis environment. No outside diagnostic-vendor results are used.

Computational reproducibility does not establish somatic status, calibrated tumor copy number, endogenous antigen presentation or drug response. The optional full-patient dense-index cloud run remains unrun; no AWS plan upgrade or EC2/EBS resources were created.
'''
manifest={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'previous_package_sha256':hashlib.sha256(old.read_bytes()).hexdigest(),'version':'v2','files':[]}
target=O/'Caris-computational-followups-2026-10-05-v2.zip'
def put(z,name,b):
 z.writestr(name,b);manifest['files'].append({'path':name,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 with zipfile.ZipFile(old) as prior:
  prior_manifest=json.loads(prior.read('MANIFEST.json'))
  for item in prior_manifest['files']:
   name=item['path']
   if '/' not in name:continue
   b=prior.read(name);assert hashlib.sha256(b).hexdigest()==item['sha256'];put(z,name,b)
 for name,p in sorted(entries.items()):
  b=p.read_bytes()
  if p.suffix in ('.json','.txt','.tsv','.csv','.py','.fasta','.fa'):assert (b'gene'+b'power') not in b.lower(),p
  put(z,name,b)
 put(z,'README.txt',readme.encode())
 put(z,'PACKAGE-EXCLUSIONS.json',(json.dumps(excluded,indent=2)+'\n').encode())
 z.writestr('MANIFEST.json',json.dumps(manifest,indent=2)+'\n')
with zipfile.ZipFile(target) as z:
 assert z.testzip() is None
 assert len(z.namelist())==len(set(z.namelist()))
 for f in manifest['files']:assert hashlib.sha256(z.read(f['path'])).hexdigest()==f['sha256']
D=L/'oct5-final-handoff-v2';D.mkdir(exist_ok=True)
top=[p for p in entries.values() if p.parent==O or p.parent==R/'output/pdf']+[target]
records=[]
for p in top:
 dest=D/p.name;shutil.copy2(p,dest);h=hashlib.sha256(p.read_bytes()).hexdigest();assert h==hashlib.sha256(dest.read_bytes()).hexdigest()
 records.append({'name':p.name,'bytes':p.stat().st_size,'sha256':h,'path':str(p),'resident':str(dest)})
receipt={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':records,'archive_member_count':len(manifest['files'])+1,'archive_members_all_SHA256_verified':True,'original_package_preserved':True,'resident_directory':str(D)}
for p in (O/'Caris-final-package-verification-2026-10-05-v2.json',D/'Caris-final-package-verification-2026-10-05-v2.json'):p.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt,indent=2))
