"""Execute only the reviewed prespecified public-read matrix; never alter filters to pass."""
from pathlib import Path
import subprocess, json, hashlib, datetime, time, os, signal, shutil, importlib.util, argparse, gzip
import pysam
P=Path(__file__).resolve().parent
C=Path('/Users/burhanazeem/.local/share/codex/caris-analysis');R=C/'genome-fusion-reference';T=C/'genome-fusion-tools'
spec=importlib.util.spec_from_file_location('original_align',P.parent/'oct4-followup/genome-fusion/align.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def run_processes(procs,status,allow_nonzero=False):
    peak=0;t=time.monotonic()
    while any(p.poll() is None for p in procs):
        free=shutil.disk_usage(P).free;rss=0
        for p in procs:
            if p.poll() is None:
                try:rss+=int(subprocess.check_output(['ps','-o','rss=','-p',str(p.pid)],text=True).strip())*1024
                except subprocess.CalledProcessError:pass
        peak=max(peak,rss)
        if free<3*1024**3 or rss>12*1024**3:
            for p in procs:
                if p.poll() is None:os.killpg(p.pid,signal.SIGTERM)
            raise RuntimeError(f'Resource guard disk={free}, RSS={rss}')
        if not allow_nonzero and any(p.poll() not in [None,0] for p in procs):raise RuntimeError('Subprocess failed '+str([p.poll() for p in procs]))
        time.sleep(.5)
    codes=[p.wait() for p in procs]
    if not allow_nonzero:assert all(c==0 for c in codes),codes
    return {'elapsed_seconds':round(time.monotonic()-t,3),'peak_observed_process_RSS_bytes':peak,'exit_codes':codes}

def run_case(d,sparse,whole=False,no_prior=False):
    label=('fullD8' if whole else f'compactD{sparse}')+('__no_known_recovery' if no_prior else '')+'__'+d['id']
    out=P/'runs'/label
    if (out/'completion.json').exists():
        done=json.loads((out/'completion.json').read_text());assert done['status']=='complete';return done
    assert not out.exists(),'Unfinished directory requires inspection: '+str(out)
    out.mkdir(parents=True)
    index=R/'STAR_G37_primary_D8_SA12' if whole else R/f'compact_BCR_ABL1/index_D{sparse}_SA12'
    reference=R/'GRCh38.primary_assembly.genome.fa' if whole else R/'compact_BCR_ABL1/assembly.fa'
    gtf=R/'gencode.v37.annotation.gtf' if whole else R/'compact_BCR_ABL1/annotation.gtf'
    if whole and not gtf.exists():
        # Use the actual uncompressed file recorded in the completed original full run.
        original=json.loads((P.parent/'oct4-followup/genome-fusion/patient-D8/run.json').read_text())['Arriba_command']
        gtf=Path(original[original.index('-g')+1])
    cmd=[str(m.STAR),'--runThreadN','2','--genomeDir',str(index),'--readFilesIn',*[x['path'] for x in d['files']],'--outFileNamePrefix',str(out/'STAR.')]+m.COMMON
    acmd=[str(m.ARRIBA),'-x',str(out/'alignments.bam'),'-o',str(out/'fusions.tsv'),'-O',str(out/'fusions.discarded.tsv'),'-a',str(reference),'-g',str(gtf),'-b',str(m.DB/'blacklist_hg38_GRCh38_v2.5.1.tsv.gz')]
    if not no_prior:acmd+=['-k',str(m.DB/'known_fusions_hg38_GRCh38_v2.5.1.tsv.gz')]
    acmd+=['-t',str(m.DB/'known_fusions_hg38_GRCh38_v2.5.1.tsv.gz'),'-p',str(m.DB/'protein_domains_hg38_GRCh38_v2.5.1.gff3')]
    metadata={'id':label,'dataset':d['id'],'started_utc':now(),'status':'running','STAR_command':cmd,'Arriba_command':acmd,'known_fusion_recovery':not no_prior,'whole_reference':whole,'sparseD':sparse,'manifest_sha256':sha(P/'design-manifest.json'),'runner_sha256':sha(Path(__file__)),'sequence_or_quality_modification':'Only as frozen in dataset manifest; runner makes none.'}
    (out/'status.json').write_text(json.dumps(metadata,indent=2))
    with (out/'alignments.bam').open('wb') as bam,(out/'STAR.stderr.log').open('w') as err:
        star=subprocess.Popen(cmd,stdout=bam,stderr=err,start_new_session=True)
        metadata['STAR_resources']=run_processes([star],metadata)
    final={}
    for line in (out/'STAR.Log.final.out').read_text().splitlines():
        if '|' in line:k,v=line.split('|',1);final[k.strip()]=v.strip()
    assert int(final['Number of input reads'])==d['input_pairs']
    metadata['STAR_final']=final
    with (out/'Arriba.log.gz').open('wb') as gz:
        a=subprocess.Popen(acmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,start_new_session=True)
        z=subprocess.Popen(['/usr/bin/gzip','-1c'],stdin=a.stdout,stdout=gz,start_new_session=True);a.stdout.close()
        metadata['Arriba_resources']=run_processes([a,z],metadata,allow_nonzero=True)
    codes=metadata['Arriba_resources']['exit_codes']
    if codes==[1,0]:
        with gzip.open(out/'Arriba.log.gz','rt') as log:message=log.read()
        assert 'ERROR: no split reads or discordant mates found' in message, 'Unrecognized Arriba failure; preserve and inspect'
        assert int(final['Number of chimeric reads'])==0, 'No-input error despite STAR chimeric reads: inspect, do not whitelist'
        pysam.samtools.quickcheck(str(out/'alignments.bam'))
        inventory={'records':0,'supplementary':0,'SA_tag':0,'mapped_different_chromosome_mates':0};names=set()
        with pysam.AlignmentFile(str(out/'alignments.bam'),'rb') as b:
            for read in b:
                inventory['records']+=1;names.add(read.query_name)
                inventory['supplementary']+=int(read.is_supplementary)
                inventory['SA_tag']+=int(read.has_tag('SA'))
                inventory['mapped_different_chromosome_mates']+=int(not read.is_unmapped and not read.mate_is_unmapped and read.reference_id!=read.next_reference_id)
        assert names==set(d['names'])
        assert inventory['supplementary']==inventory['SA_tag']==inventory['mapped_different_chromosome_mates']==0, inventory
        junction_rows=[line for line in (out/'STAR.Chimeric.out.junction').read_text().splitlines() if line and not line.startswith('#') and not line.startswith('chr_donorA')]
        assert not junction_rows
        inventory['query_names']=len(names);inventory['chimeric_junction_rows']=0;inventory['BAM_quickcheck']='passed'
        metadata['no_chimeric_input_validation']=inventory
        metadata['caller_outcome']='no_chimeric_input; stock Arriba exit1, no call table generated'
        metadata['caller_error_message']=message.split('ERROR:',1)[1].strip()
        assert not (out/'fusions.tsv').exists() and not (out/'fusions.discarded.tsv').exists()
    else:
        assert codes==[0,0],codes
        metadata['caller_outcome']='completed_with_call_tables'
        assert (out/'fusions.tsv').exists() and (out/'fusions.discarded.tsv').exists()
    metadata.update(status='complete',finished_utc=now(),outputs={f.name:{'sha256':sha(f),'bytes':f.stat().st_size} for f in out.iterdir() if f.is_file() and f.name not in ['status.json','completion.json']})
    (out/'status.json').write_text(json.dumps(metadata,indent=2));(out/'completion.json').write_text(json.dumps(metadata,indent=2))
    print(json.dumps({'completed':label,'input_pairs':d['input_pairs'],'STAR_seconds':metadata['STAR_resources']['elapsed_seconds'],'Arriba_seconds':metadata['Arriba_resources']['elapsed_seconds']}),flush=True)
    return metadata

if __name__=='__main__':
    design=json.loads((P/'design-manifest.json').read_text());freeze=json.loads((P/'design-freeze.json').read_text())
    assert sha(P/'design-manifest.json')==freeze['manifest_sha256']
    review=json.loads((P/'independent-design-review.json').read_text())
    assert review['status']=='approved' and review['manifest_sha256']==freeze['manifest_sha256'],'Independent review must approve this exact design'
    for d in design['datasets']:
        for x in d['files']:assert sha(x['path'])==x['sha256']
    datasets={d['id']:d for d in design['datasets']}
    (P/'execution-start.json').write_text(json.dumps({'started_utc':now(),'design_sha256':freeze['manifest_sha256'],'runner_sha256':sha(Path(__file__)),'review_sha256':sha(P/'independent-design-review.json')},indent=2))
    ca=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())]);done=[]
    try:
        for case in design['compact_matrix']:done.append(run_case(datasets[case['dataset']],case['sparseD']))
        for name in design['full_reference_D8_transfer_checks']:done.append(run_case(datasets[name],8,whole=True))
        x=design['no_prior_comparison'];done.append(run_case(datasets[x['dataset']],x['sparseD'],no_prior=True))
        (P/'execution-complete.json').write_text(json.dumps({'finished_utc':now(),'completed_runs':len(done),'design_sha256':freeze['manifest_sha256'],'all_exit_codes_zero':all(x['Arriba_resources']['exit_codes']==[0,0] for x in done),'no_chimeric_input_exit1_cases':[x['id'] for x in done if x['Arriba_resources']['exit_codes']==[1,0]],'run_ids':[x['id'] for x in done]},indent=2))
    finally:ca.terminate()
