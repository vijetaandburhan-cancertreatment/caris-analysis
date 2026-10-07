"""Bounded local-only Pizzly research workflow; public setup must finish first."""
import gzip, json, pathlib, subprocess, sys
import pysam

root = pathlib.Path(__file__).resolve().parent
source = pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
stage = sys.argv[1]
kallisto = root / 'public/kallisto-v0.44.0/kallisto/kallisto'
pizzly = root / 'public/pizzly-v0.37.3/pizzly'
fasta = root / 'public/gencode.v37.transcripts.pizzly-headers.fa.gz'
gtf = root.parent / 'genomics/gencode.v37.annotation.gtf.gz'
index = root / 'public/gencode.v37.k31.kallisto044.idx'
raw = [source / f'RNA_TN26-279853_S25.R{i}.fastq.gz' for i in [1,2]]

def run(name, cmd):
    subprocess.run([sys.executable, str(root/'run_guarded.py'), name, *[str(x) for x in cmd]], check=True)

if stage == 'index':
    assert (root/'public-input-manifest.json').exists(), 'Public reference acquisition incomplete'
    run('gencode37-index-attempt2', [kallisto, 'index', '-k', '31', '-i', index, fasta])
elif stage == 'public-full-reference':
    out=root/'public-full-reference';out.mkdir(exist_ok=True)
    reads=[root/f'public/pizzly-v0.37.3/test/reads_{i}.fastq.gz' for i in [1,2]]
    run('public-full-reference-quant',[kallisto,'quant','-t','2','--fusion','-i',index,'-o',out/'kallisto',*reads])
    run('public-full-reference-pizzly',[pizzly,'-k','31','--gtf',gtf,'--cache',root/'public/gencode37.pizzly.cache','--align-score','2','--insert-size','400','--fasta',fasta,'--output',out/'pizzly',out/'kallisto/fusion.txt'])
    expected={tuple(sorted(x)) for x in [('EML4','ALK'),('CD74','ROS1'),('HOOK3','RET'),('TMPRSS2','ETV1'),('EWSR1','FLI1'),('ETV6','NTRK3'),('BRD4','NUTM1'),('EWSR1','ATF1'),('AKAP9','BRAF')]}
    observed={tuple(sorted((x['geneA']['name'],x['geneB']['name']))) for x in json.loads((out/'pizzly.json').read_text())['genes']}
    result={'expected_public_example_pairs':sorted(expected),'observed':sorted(observed),'recovered':sorted(expected&observed),'missing':sorted(expected-observed),'interpretation':'Checks reference/header/annotation interoperability on supplied public positive examples; not specimen sensitivity validation'}
    (out/'control-summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    assert len(expected&observed)>=5,'Many expected positive controls absent: stop and investigate'
elif stage == 'smoke':
    out = root / 'patient-smoke'
    out.mkdir(exist_ok=True)
    paired = [out/'RNA.sample.R1.fastq.gz', out/'RNA.sample.R2.fastq.gz']
    n = 0
    with pysam.FastxFile(str(raw[0])) as r1, pysam.FastxFile(str(raw[1])) as r2, gzip.open(paired[0],'wt',compresslevel=1) as w1, gzip.open(paired[1],'wt',compresslevel=1) as w2:
        for a,b in zip(r1,r2):
            assert a.name == b.name
            for r,w in [(a,w1),(b,w2)]:
                header = r.name + (' '+r.comment if r.comment else '')
                w.write('@'+header+'\n'+r.sequence+'\n+\n'+r.quality+'\n')
            n += 1
            if n == 100000:
                break
    (out/'sampling.json').write_text(json.dumps({'paired_reads':n,'selection':'First100000 input pairs for smoke test only, not sensitivity estimate','sources':list(map(str,raw))},indent=2)+'\n')
    run('patient-smoke-quant', [kallisto,'quant','-t','2','--fusion','-i',index,'-o',out/'kallisto',*paired])
    info=json.loads((out/'kallisto/run_info.json').read_text())
    assert info['n_processed']==n
    assert (out/'kallisto/fusion.txt').exists()
    run('patient-smoke-pizzly', [pizzly,'-k','31','--gtf',gtf,'--cache',root/'public/gencode37.pizzly.cache','--align-score','2','--insert-size','400','--fasta',fasta,'--output',out/'pizzly',out/'kallisto/fusion.txt'])
    print(json.dumps({'passed':True,'run_info':info,'pizzly_pairs':len(json.loads((out/'pizzly.json').read_text())['genes'])},indent=2))
elif stage == 'full-quant':
    out=root/'full';out.mkdir(exist_ok=True)
    assert (root/'patient-smoke/pizzly.json').exists()
    run('full-quant',[kallisto,'quant','-t','2','--fusion','-i',index,'-o',out/'kallisto',*raw])
    info=json.loads((out/'kallisto/run_info.json').read_text())
    assert info['n_processed']==23209264,(info['n_processed'],'Full pairing count mismatch')
    assert (out/'kallisto/fusion.txt').exists()
    print(json.dumps(info,indent=2))
elif stage == 'full-pizzly':
    insert=int(sys.argv[2]); out=root/'full'
    run(f'full-pizzly-insert{insert}',[pizzly,'-k','31','--gtf',gtf,'--cache',root/'public/gencode37.pizzly.cache','--align-score','2','--insert-size',str(insert),'--fasta',fasta,'--output',out/f'pizzly.insert{insert}',out/'kallisto/fusion.txt'])
else:
    raise SystemExit('Stage must be index, public-full-reference, smoke, full-quant or full-pizzly INSERT_SIZE')
