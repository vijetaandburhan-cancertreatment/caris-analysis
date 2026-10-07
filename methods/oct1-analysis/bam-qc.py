import json
import pathlib
import time
import pysam

ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'outputs/caris-raw-data/TN26-279853'
OUT = ROOT / 'work/oct1-analysis'
results = {'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'pysam_version': pysam.__version__, 'samtools_version': pysam.__samtools_version__, 'max_threads': 2, 'files': []}
for kind in ('RNA', 'DNA'):
    path = SOURCE / f'{kind}_TN26-279853.bam'
    index = OUT / f'{path.name}.bai'
    started = time.time()
    print(f'{kind}: indexing {path.stat().st_size} bytes', flush=True)
    pysam.index('-@', '2', '-o', str(index), str(path))
    print(f'{kind}: flagstat', flush=True)
    stats = json.loads(pysam.flagstat('-@', '2', '-O', 'json', str(path)))
    with pysam.AlignmentFile(str(path), 'rb', index_filename=str(index), threads=2) as bam:
        idxstats = [{'contig': item.contig, 'mapped': item.mapped, 'unmapped': item.unmapped, 'total': item.total} for item in bam.get_index_statistics()]
        header = bam.header.to_dict()
        item = {'kind': kind, 'path': str(path), 'bytes': path.stat().st_size, 'index': str(index), 'index_bytes': index.stat().st_size, 'header': header, 'flagstat': stats, 'idxstats': idxstats, 'elapsed_seconds': round(time.time()-started,1)}
        results['files'].append(item)
    (OUT / 'bam-qc.json').write_text(json.dumps(results, indent=2) + '\n')
    print(f'{kind}: complete {item["elapsed_seconds"]} sec; {json.dumps(stats)}', flush=True)
print('Complete', flush=True)
