"""Small, local BAM->paired FASTQ extraction checks; no full patient extraction.

RNA HLA-A locus plus1kb is only a mechanics test, not a valid complete HLA panel
input. Unmapped/all-chr6 selection must be separately done for real inference.
"""
from pathlib import Path
import collections,gzip,hashlib,json,re,time
import pysam
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
ARC=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/hla-tools/arcasHLA')
SRC=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
TMP=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/hla-tools/smoke')
TMP.mkdir(exist_ok=True)
hla=None
for line in gzip.open(ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz','rt'):
    if '\tgene\t' in line and 'gene_name "HLA-A"' in line:
        a=line.rstrip().split('\t');hla=(a[0],int(a[3])-1,int(a[4]));break
assert hla is not None
result=[]
for label,path,index,region in [
    ('public_test',ARC/'test/test.bam',None,None),
    ('patient_HLA_A_mechanics_only',SRC/'RNA_TN26-279853.bam',ROOT/'work/oct1-analysis/RNA_TN26-279853.bam.bai',(hla[0],hla[1]-1000,hla[2]+1000))]:
    seen=collections.defaultdict(dict);counts=collections.Counter();fingerprints={};duplicate_primary_keys=0
    kwargs={'index_filename':str(index)} if index else {}
    with pysam.AlignmentFile(str(path),'rb',**kwargs) as b:
        it=b.fetch(*region) if region else b.fetch(until_eof=True)
        for r in it:
            counts['input_alignments']+=1
            if r.flag&(4|256|2048) or not r.is_proper_pair:counts['filtered_not_primary_proper_pair']+=1;continue
            if r.query_sequence is None:counts['missing_sequence']+=1;continue
            if r.query_qualities is None and label!='public_test':counts['missing_patient_quality_excluded']+=1;continue
            mate=1 if r.is_read1 else 2 if r.is_read2 else 0
            if not mate:counts['no_mate_number']+=1;continue
            k=(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name)
            if mate in seen[k]:duplicate_primary_keys+=1;continue
            seq=r.get_forward_sequence();qual=r.get_forward_qualities() if r.query_qualities is not None else None
            if qual is None:
                # Upstream public test BAM intentionally has no qualities.
                # Kallisto uses sequence, not these placeholders. Never used
                # for patient extraction or quantitative base-quality claims.
                qual=[40]*len(seq);counts['public_only_placeholder_quality_reads']+=1
            assert len(seq)==len(qual)
            seen[k][mate]=(seq,''.join(chr(q+33) for q in qual));counts['eligible_reads']+=1
    assert duplicate_primary_keys==0
    files=[TMP/(label+'.R1.fq.gz'),TMP/(label+'.R2.fq.gz')]
    digest=[hashlib.sha256(),hashlib.sha256()];n=0
    with gzip.open(files[0],'wt') as f1,gzip.open(files[1],'wt') as f2:
        for (rg,name),mates in sorted(seen.items()):
            if set(mates)!={1,2}:counts['orphan_query_names_not_written']+=1;continue
            # RG is included so a collision across read groups cannot merge pairs.
            qname=name if not rg else rg+'|'+name
            for m,f in [(1,f1),(2,f2)]:
                seq,qual=mates[m];s='@'+qname+'\n'+seq+'\n+\n'+qual+'\n';f.write(s);digest[m-1].update(s.encode())
            n+=1
    recount=0;rd=[hashlib.sha256(),hashlib.sha256()]
    with pysam.FastxFile(str(files[0])) as f1,pysam.FastxFile(str(files[1])) as f2:
        i1=iter(f1);i2=iter(f2)
        while True:
            a=next(i1,None);b=next(i2,None)
            if a is None or b is None:assert a is None and b is None;break
            assert a.name==b.name
            for j,z in enumerate((a,b)):
                assert len(z.sequence)==len(z.quality)
                rd[j].update(('@'+z.name+'\n'+z.sequence+'\n+\n'+z.quality+'\n').encode())
            recount+=1
    assert n>0 and recount==n and [d.hexdigest() for d in digest]==[d.hexdigest() for d in rd]
    result.append({'label':label,'source':str(path),'region0':region,'counts':dict(counts),'complete_pairs_written_and_roundtrip_verified':n,'duplicate_primary_keys':duplicate_primary_keys,'paired_names_and_sequence_quality_exact_roundtrip':True,'FASTQ_content_sha256':[d.hexdigest() for d in rd],'gzip_bytes':[p.stat().st_size for p in files],'outputs':[str(p) for p in files]})
(OUT/'extraction-smoke.json').write_text(json.dumps({'status':'passed','method':'Primary proper-paired alignments, no MAPQ threshold/PCR deduplication. Reverse-strand reads restored to original sequence/quality orientation. Only complete RG+query-name pairs written. Orphans counted separately. No read names or sequences logged.','limitations':'Patient HLA-A locus subset is a mechanics control; never use it to make whole-panel HLA genotypes. Public test BAM lacks quality strings: public-test-only Q40 placeholder strings supplied for sequence-only kallisto, never patient data. Public test BAM extraction is proper-paired all-read subset.','results':result},indent=2)+'\n')
print(json.dumps(result,indent=2))
