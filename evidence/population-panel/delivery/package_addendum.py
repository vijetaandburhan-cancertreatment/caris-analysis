from pathlib import Path
import hashlib,json,zipfile,shutil,datetime
D=Path(__file__).resolve().parent;R=D.parent
W=Path('/Users/burhanazeem/Documents/Codex/2026-09-05/finances-plugin-finances-openai-curated-remote-3')
O=W/'outputs/caris-followup'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def js(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
assert json.loads((D/'final-content-approval.json').read_text())['status']=='PASS'
assert json.loads((D/'editorial.json').read_text())['ready_for_final_render']
assert json.loads((R/'independent-count-audit/DNA-audit.json').read_text())['status']=='PASS'
assert json.loads((R/'independent-count-audit/RNA-audit.json').read_text())['status']=='PASS'
assert json.loads((R/'independent-count-audit/independent-summary-audit.json').read_text())['status']=='PASS'
prior=R.parent/'oct5-final-handoff-v2'
expected={'Caris-DNA-RNA-research-handoff-TN26-279853-2026-10-05-v2.pdf':'d1ca641a0235b7b6b09846c7e2b870dea36dbd6e5cf17ca0f4709c3eec784c9d','Caris-followup-findings-2026-10-05-v2.pdf':'0a851414a3235d829ee4a8c7e423feed307c85a5ffdf865ee61011adf3267c03','Caris-computational-followups-2026-10-05-v2.zip':'3d344e7f4722353217d3a13cbd74fa5ec3cd6667132d1520232e972fb9d975d8'}
oldchecks={}
for n,h in expected.items():
 assert sha(prior/n)==h,n;oldchecks[n]=h
pdf=D/'Caris-DNA-RNA-population-panel-addendum-2026-10-05.pdf';assert pdf.exists()
archive=D/'Caris-population-panel-evidence-2026-10-05.zip'
excluded=[];included=[];skipnames={'package-manifest.json','package-verification.json','package-checksums.sha256',archive.name}
for p in sorted(R.rglob('*')):
 if not p.is_file():continue
 rel=p.relative_to(R)
 if '__pycache__' in rel.parts or any('preview' in s for s in rel.parts) or 'final-render' in rel.parts or p.name in skipnames or p.suffix=='.log':continue
 if p.suffix=='.zip':continue
 rec={'path':str(rel),'bytes':p.stat().st_size,'sha256':sha(p)}
 if p.name.endswith('.pileup.gz') and rel.parts[0] in ('DNA','RNA'):excluded.append(rec)
 else:included.append(rec)
manifest={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Completed population-panel technical addendum; prior v2 milestone unchanged','included_files':included,'excluded_large_raw_pileups_preserved_resident':excluded,'resident_root':str(R),'unchanged_v2':oldchecks,'sources':'Original Caris DNA/RNA and public reference; no external patient-data uploads or clinical outreach in this follow-up.'}
js(D/'package-manifest.json',manifest)
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for r in included:z.write(R/r['path'],arcname=r['path'])
 z.write(D/'package-manifest.json',arcname='delivery/package-manifest.json')
with zipfile.ZipFile(archive) as z:
 assert z.testzip() is None
 for r in included:assert hashlib.sha256(z.read(r['path'])).hexdigest()==r['sha256'],r['path']
 assert hashlib.sha256(z.read('delivery/package-manifest.json')).hexdigest()==sha(D/'package-manifest.json')
files=[pdf,D/(pdf.stem+'.txt'),archive,D/'package-manifest.json',D/'Caris-population-panel-morning-update-2026-10-05.txt']
copychecks=[]
for p in files:
 dest=O/p.name;assert not dest.exists(),f'Will not overwrite {dest}'
 shutil.copy2(p,dest);assert sha(dest)==sha(p)
 copychecks.append({'path':str(dest),'resident_path':str(p),'bytes':p.stat().st_size,'sha256':sha(p)})
verification={'status':'PASS','archive_members':len(included)+1,'member_hashes_verified':len(included)+1,'excluded_raw_pileups_preserved':len(excluded),'outputs':copychecks,'unchanged_v2':oldchecks,'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
js(D/'package-verification.json',verification);shutil.copy2(D/'package-verification.json',O/'Caris-population-panel-package-verification-2026-10-05.json')
(D/'package-checksums.sha256').write_text('\n'.join(r['sha256']+'  '+Path(r['resident_path']).name for r in copychecks)+'\n')
print(json.dumps(verification,indent=2))
