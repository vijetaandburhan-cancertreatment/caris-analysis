"""Restore public models from pinned manifests; never reads patient input.
Saves next to this script; refuses checksum-mismatched existing files.
"""
import pathlib,json,hashlib,urllib.request
P=pathlib.Path(__file__).resolve().parent
for manifest in ['CapHLA-source-manifest.json','CapHLA-v1-source-manifest.json']:
 for r in json.loads((P/manifest).read_text())['files']:
  rel=pathlib.Path(r['path']).relative_to('work/oct4-followup/hla-II'); f=P/rel; f.parent.mkdir(parents=True,exist_ok=True)
  if f.exists():
   assert hashlib.sha256(f.read_bytes()).hexdigest()==r['sha256'],str(f)
   continue
  b=urllib.request.urlopen(r['url']).read();assert len(b)==r['bytes'];assert hashlib.sha256(b).hexdigest()==r['sha256'];f.write_bytes(b)
