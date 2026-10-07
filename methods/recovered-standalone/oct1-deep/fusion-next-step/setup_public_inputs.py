"""Download only pinned public software and GENCODE transcript reference, never patient data."""
import datetime, gzip, hashlib, json, pathlib, tarfile, urllib.request

ROOT = pathlib.Path(__file__).resolve().parent
PUBLIC = ROOT / 'public'
PUBLIC.mkdir(exist_ok=True)
URLS = {
    'kallisto_mac_m1-v0.51.1.tar.gz': 'https://github.com/pachterlab/kallisto/releases/download/v0.51.1/kallisto_mac_m1-v0.51.1.tar.gz',
    'pizzly_mac-v0.37.3.tar.gz': 'https://github.com/pmelsted/pizzly/releases/download/v0.37.3/pizzly_mac.tar.gz',
    'gencode.v37.transcripts.fa.gz': 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz',
}
manifest = []
for name, url in URLS.items():
    p = PUBLIC / name
    if not p.exists():
        with urllib.request.urlopen(url, timeout=60) as r, open(p, 'wb') as f:
            n = 0
            while b := r.read(1024 * 1024):
                n += len(b)
                if n > 100_000_000:
                    raise RuntimeError('Per-file public download cap exceeded')
                f.write(b)
    b = p.read_bytes()
    item = {'url': url, 'path': str(p), 'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest()}
    if name.startswith(('kallisto_', 'pizzly_')):
        with tarfile.open(p, 'r:gz') as tf:
            item['archive_members'] = [m.name for m in tf.getmembers()]
            dest = PUBLIC / ('kallisto-v0.51.1' if name.startswith('kallisto') else 'pizzly-v0.37.3')
            dest.mkdir(exist_ok=True)
            for m in tf.getmembers():
                q = (dest / m.name).resolve()
                if not q.is_relative_to(dest.resolve()) or m.issym() or m.islnk():
                    raise RuntimeError('Unsafe archive member')
            tf.extractall(dest, filter='data')
    manifest.append(item)

# GENCODE headers use pipes; Pizzly's official documentation requires replacing these
# in the same FASTA supplied both to kallisto indexing and Pizzly annotation.
fixed = PUBLIC / 'gencode.v37.transcripts.pizzly-headers.fa.gz'
seq_count = seq_bases = 0
with gzip.open(PUBLIC / 'gencode.v37.transcripts.fa.gz', 'rb') as src, open(fixed, 'wb') as raw:
    with gzip.GzipFile(filename='', mode='wb', fileobj=raw, compresslevel=6, mtime=0) as dst:
        for line in src:
            if line.startswith(b'>'):
                line = line.replace(b'|', b' ')
                seq_count += 1
            else:
                seq_bases += len(line.strip())
            dst.write(line)
manifest.append({'path': str(fixed), 'bytes': fixed.stat().st_size,
                 'sha256': hashlib.sha256(fixed.read_bytes()).hexdigest(),
                 'transform': 'Replace pipe delimiters with spaces in FASTA headers only; sequences unmodified',
                 'sequences': seq_count, 'sequence_bases': seq_bases})
(ROOT / 'public-input-manifest.json').write_text(json.dumps({'time_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'files': manifest}, indent=2) + '\n')
print(json.dumps(manifest, indent=2))
