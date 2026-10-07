"""Use upstream Pizzly's public test data to test pinned-tool interoperability first."""
import json, os, pathlib, subprocess, sys

root = pathlib.Path(__file__).resolve().parent
test = root / 'public/pizzly-v0.37.3/test'
out = root / os.environ.get('SMOKE_OUTPUT_NAME', 'public-smoke')
out.mkdir(exist_ok=True)
kallisto = pathlib.Path(os.environ.get('KALLISTO_PATH', str(root / 'public/kallisto-v0.51.1/kallisto/kallisto')))
pizzly = root / 'public/pizzly-v0.37.3/pizzly'
commands = [
 [out.name+'-index', str(kallisto), 'index', '-k', '31', '-i', str(out/'transcripts.idx'), str(test/'transcripts.fasta.gz')],
 [out.name+'-fusion', str(kallisto), 'quant', '-t', '2', '--fusion', '-i', str(out/'transcripts.idx'), '-o', str(out/'kallisto'), str(test/'reads_1.fastq.gz'), str(test/'reads_2.fastq.gz')],
 [out.name+'-pizzly', str(pizzly), '-k', '31', '--gtf', str(test/'transcripts.gtf.gz'), '--cache', str(out/'annotation.cache'), '--align-score', '2', '--insert-size', '400', '--fasta', str(test/'transcripts.fasta.gz'), '--output', str(out/'pizzly'), str(out/'kallisto/fusion.txt')],
]
for cmd in commands:
    subprocess.run([sys.executable, str(root/'run_guarded.py'), *cmd], check=True)
summary = {'purpose':'Upstream public test for binary, HDF5 ABI, gzipped annotation and kallisto/Pizzly output compatibility. Not a clinical analytical validation.', 'kallisto_run_info': json.loads((out/'kallisto/run_info.json').read_text()), 'fusion_txt_bytes':(out/'kallisto/fusion.txt').stat().st_size, 'pizzly_output':json.loads((out/'pizzly.json').read_text())}
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2)[:6000])
