"""Read-only, bounded-memory structural scan; no sequence data is written."""
import pathlib,json,time,collections,itertools,traceback,sys,os
import pysam
BASE=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
OUT=pathlib.Path(__file__).resolve().parent
CHECK=OUT/'fastq-structural-qc.checkpoint.json'
RESULT=OUT/'fastq-structural-qc.json'
BAM=json.loads((OUT/'bam-qc.json').read_text())
EXPECTED={r['kind']:r['flagstat']['QC-passed reads']['read1'] for r in BAM['files']}
allow_slow='--continue-slow' in sys.argv
results={'pysam_version':pysam.__version__,'method':'Sequential paired streaming using pysam.FastxFile, every name pair and sequence/quality length checked; no extracted FASTQ written. First1000 pairs and then one pair per10000 sampled for quality/length histograms. No raw sequence or read names recorded. One process, default single-thread decompression.','started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'status':'running','pairs':[]}

def save(status, current=None):
    d={**results,'status':status,'updated_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
    if current:d['current']=current
    tmp=CHECK.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2)+'\n');tmp.replace(CHECK)

def norm(name):
    return name[:-2] if name.endswith(('/1','/2')) else name

try:
    for kind in ('DNA','RNA'):
        t0=time.monotonic(); n=0; mismatches=0; length_bad=0; missing_quality=0; empty=0
        minlen=[None,None];maxlen=[0,0];basecounts=[0,0]
        qs=[collections.Counter(),collections.Counter()];ls=[collections.Counter(),collections.Counter()]
        sample_n=0;errors=[]
        paths=[BASE/f'{kind}_TN26-279853_S25.R{i}.fastq.gz' for i in (1,2)]
        size=[p.stat().st_size for p in paths]
        print(json.dumps({'event':'start','kind':kind,'expected_pairs':EXPECTED[kind],'compressed_bytes':size}),flush=True)
        f1=pysam.FastxFile(str(paths[0]),persist=False);f2=pysam.FastxFile(str(paths[1]),persist=False)
        exhausted=False
        try:
            for a,b in itertools.zip_longest(f1,f2):
                n+=1
                if a is None or b is None:
                    errors.append({'pair_ordinal':n,'error':'unequal_mate_record_counts','missing_mate':1 if a is None else 2});
                    raise ValueError('FASTQ mate record counts unequal')
                if norm(a.name)!=norm(b.name):
                    mismatches+=1
                    if len(errors)<10:errors.append({'pair_ordinal':n,'error':'mate_name_mismatch'})
                sample=n<=1000 or n%10000==0
                if sample:sample_n+=1
                for i,r in enumerate((a,b)):
                    seq=r.sequence;qual=r.quality;l=len(seq)
                    basecounts[i]+=l;minlen[i]=l if minlen[i] is None else min(minlen[i],l);maxlen[i]=max(maxlen[i],l)
                    if l==0:empty+=1
                    if qual is None:
                        missing_quality+=1
                        if len(errors)<10:errors.append({'pair_ordinal':n,'mate':i+1,'error':'missing_quality'})
                    elif len(qual)!=l:
                        length_bad+=1
                        if len(errors)<10:errors.append({'pair_ordinal':n,'mate':i+1,'error':'seq_quality_length_mismatch'})
                    if sample:
                        ls[i][l]+=1
                        if qual is not None:qs[i].update(qual)
                if n%1000000==0:
                    elapsed=time.monotonic()-t0
                    current={'kind':kind,'pairs_scanned':n,'expected_BAM_primary_pairs':EXPECTED[kind],'elapsed_seconds':round(elapsed,2),'estimated_total_seconds_from_rate':round(elapsed*EXPECTED[kind]/n,2),'pair_name_mismatches':mismatches,'seq_quality_length_mismatches':length_bad,'missing_quality':missing_quality}
                    save('running',current);print(json.dumps({'event':'progress',**current}),flush=True)
                    if n==1000000 and current['estimated_total_seconds_from_rate']>900 and not allow_slow:
                        save('paused_after_benchmark_over15minutes',current)
                        print(json.dumps({'event':'benchmark_pause',**current}),flush=True)
                        sys.exit(75)
            exhausted=True
        finally:
            f1.close();f2.close()
        elapsed=time.monotonic()-t0
        qsum=[]
        for counts in qs:
            total=sum(counts.values());hist={str(ord(k)-33):v for k,v in sorted(counts.items())}
            qsum.append({'sampled_quality_characters':total,'phred33_histogram':hist,'fraction_Q20_or_above':sum(v for k,v in counts.items() if ord(k)-33>=20)/total if total else None,'fraction_Q30_or_above':sum(v for k,v in counts.items() if ord(k)-33>=30)/total if total else None})
        r={'kind':kind,'source_paths':[str(p) for p in paths],'source_bytes':size,'records_per_mate':n,'normal_EOF_both':exhausted,'matched_BAM_primary_pair_count':n==EXPECTED[kind],'BAM_primary_pairs':EXPECTED[kind],'pair_name_mismatches':mismatches,'seq_quality_length_mismatches':length_bad,'missing_quality':missing_quality,'zero_length_records':empty,'minimum_read_lengths':minlen,'maximum_read_lengths':maxlen,'total_bases_per_mate':basecounts,'sampled_pairs':sample_n,'sampled_length_histograms':ls,'sampled_quality':qsum,'error_examples_without_sequence':errors,'elapsed_seconds':round(elapsed,2),'status':'passed' if exhausted and mismatches==length_bad==missing_quality==empty==0 and n==EXPECTED[kind] else 'needs_review','gzip_integrity_note':'Both compressed FASTQ streams consumed to normal parser EOF with no decompression/parser exception; original compressed bytes separately source-checksum verified.'}
        results['pairs'].append(r);save('running');RESULT.write_text(json.dumps(results,indent=2)+'\n');print(json.dumps({'event':'complete','kind':kind,'records_per_mate':n,'elapsed_seconds':r['elapsed_seconds'],'status':r['status']}),flush=True)
    results['status']='passed' if all(x['status']=='passed' for x in results['pairs']) else 'needs_review';results['finished_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime());RESULT.write_text(json.dumps(results,indent=2)+'\n');save(results['status']);print(json.dumps({'event':'all_complete','status':results['status']}),flush=True)
except Exception as e:
    results['status']='error';results['exception_type']=type(e).__name__;results['exception_message']=str(e);RESULT.write_text(json.dumps(results,indent=2)+'\n');save('error',{'kind':kind,'pairs_scanned':n});print(json.dumps({'event':'error','kind':kind,'pairs_scanned':n,'error_type':type(e).__name__,'message':str(e)}),flush=True);raise
