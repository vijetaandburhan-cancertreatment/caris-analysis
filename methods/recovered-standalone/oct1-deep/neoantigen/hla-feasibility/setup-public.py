"""Fetch immutable official HLA sources/reference, no patient network transfer."""
from pathlib import Path
import hashlib,json,urllib.request,time
OUT=Path(__file__).resolve().parent
DEST=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/hla-tools')
ARCAS='9fa54a212d134b0d9894d1fc19ec1bdc6f62eb55'
IMGT='8d77b3dd93959663d58ae5b626289d0746edd0e7'
ARC=DEST/'arcasHLA';HLA=ARC/'dat/IMGTHLA'
sources=[]
def get(url,p,expected_sha=None,maxsize=250_000_000):
    p.parent.mkdir(parents=True,exist_ok=True)
    if not p.exists():
        tmp=p.with_suffix(p.suffix+'.part');n=0
        with urllib.request.urlopen(url,timeout=60) as r,tmp.open('wb') as f:
            while True:
                buf=r.read(1024*1024)
                if not buf:break
                n+=len(buf)
                if n>maxsize:raise RuntimeError('Public file exceeds preflight cap')
                f.write(buf)
        tmp.rename(p)
    sha=hashlib.sha256(p.read_bytes()).hexdigest()
    if expected_sha:assert sha==expected_sha,(p,sha)
    sources.append({'path':str(p),'source_url':url,'bytes':p.stat().st_size,'sha256':sha})
    print(p.name,p.stat().st_size,sha,flush=True)
for name in ['LICENSE','README.md','arcasHLA','environment.yml','scripts/align.py','scripts/arcas_utilities.py','scripts/convert.py','scripts/customize.py','scripts/extract.py','scripts/genotype.py','scripts/merge.py','scripts/partial.py','scripts/quant.py','scripts/reference.py','dat/info/decoys_alts.json','dat/info/hla_freq.tsv','dat/info/parameters.json','test/test.bam','test/expected_output/test.genotype.json']:
    get(f'https://raw.githubusercontent.com/RabadanLab/arcasHLA/{ARCAS}/{name}',ARC/name,maxsize=5_000_000)
for name in ['LICENCE.md','release_version.txt','wmda/hla_nom_g.txt','wmda/hla_nom_p.txt']:
    get(f'https://raw.githubusercontent.com/ANHIG/IMGTHLA/{IMGT}/{name}',HLA/name,maxsize=5_000_000)
get(f'https://media.githubusercontent.com/media/ANHIG/IMGTHLA/{IMGT}/hla.dat',HLA/'hla.dat','939680ee68ed2b4cfb4773af6c8e0b5fc5c985de00611ba13e3fcba3197d898f')
assert sum(z['bytes'] for z in sources)<1_000_000_000
(ARC/'dat/ref').mkdir(exist_ok=True)
manifest={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'arcas_commit':ARCAS,'IMGT_commit':IMGT,'IMGT_release':'3.46.0','files':sources,'downloaded_bytes':sum(z['bytes'] for z in sources),'note':'Official catalog version explicitly represented in arcasHLA parameters. Older than current database; inference remains exploratory. No patient information in requests.'}
(OUT/'public-reference-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
