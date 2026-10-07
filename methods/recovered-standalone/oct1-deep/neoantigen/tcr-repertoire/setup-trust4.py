"""Public, pinned TRUST4 source download and bounded task-local build only."""
from pathlib import Path
import hashlib,json,subprocess,tarfile,time,urllib.request
OUT=Path(__file__).resolve().parent
BASE=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/trust4-tools');BASE.mkdir(exist_ok=True)
COMMIT='a3fedd4aa0c1ad4da815427d82e0ceb69ce9c3a0'
URL='https://codeload.github.com/liulab-dfci/TRUST4/tar.gz/'+COMMIT
archive=BASE/('TRUST4-'+COMMIT+'.tar.gz');source=BASE/('TRUST4-'+COMMIT)
assert not archive.exists() and not source.exists(), 'Preserve prior setup'
with urllib.request.urlopen(URL,timeout=45) as response,archive.open('wb') as f:
    size=0
    while True:
        b=response.read(1024*1024)
        if not b:break
        size+=len(b);assert size<100*1024**2;f.write(b)
with tarfile.open(archive) as tf:
    assert sum(x.size for x in tf.getmembers())<500*1024**2
    for x in tf.getmembers():
        assert not x.name.startswith('/') and '..' not in Path(x.name).parts
        assert not x.issym() and not x.islnk(), 'Unexpected link'
    tf.extractall(BASE,filter='data')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest={'source_url':URL,'commit':COMMIT,'archive_sha256':sha(archive),'archive_bytes':archive.stat().st_size,'source_dir':str(source),'license':'MIT (source LICENSE.txt), bundled samtools notices retained','sources':[{'file':str(p.relative_to(source)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(source.rglob('*')) if p.is_file()],'compiler':subprocess.check_output(['clang++','--version'],text=True)}
(OUT/'public-source-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
args=['make','-j','1','CXX=clang++','CXXFLAGS=-O3 -Wall -g -std=gnu++11']
start=time.time()
with (OUT/'build.log').open('w') as f:
    p=subprocess.run(args,cwd=source,stdout=f,stderr=subprocess.STDOUT,timeout=600)
build={'command':args,'returncode':p.returncode,'elapsed_seconds':time.time()-start,'directory_bytes':sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file()),'executables':[{'path':str(source/name),'sha256':sha(source/name),'bytes':(source/name).stat().st_size} for name in ['trust4','bam-extractor','fastq-extractor','annotator'] if (source/name).is_file()]}
(OUT/'build-result.json').write_text(json.dumps(build,indent=2)+'\n')
print(json.dumps(build,indent=2));assert p.returncode==0
