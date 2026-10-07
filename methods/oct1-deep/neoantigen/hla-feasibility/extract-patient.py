"""Local chr6 proper-pair + both-unmapped primary RNA selection for arcasHLA.

Preserve raw patient qualities and read orientation. Disk-backed bounded name sort;
no MAPQ cutoff, no deduplication. Never logs read names or sequences.
"""
from pathlib import Path
from collections import Counter, defaultdict
import gzip, hashlib, itertools, json, os, resource, shutil, time
import pysam

OUT=Path(__file__).resolve().parent; ROOT=OUT.parents[3]
BASE=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/hla-tools')
DEST=BASE/'patient';DEST.mkdir(exist_ok=True)
SRC=BASE.parent/'TN26-279853/RNA_TN26-279853.bam'
IDX=ROOT/'work/oct1-analysis/RNA_TN26-279853.bam.bai'
selected=DEST/'selected-primary.bam';sorted_bam=DEST/'selected-primary.namesort.bam'
counts=Counter();stages=[];started=time.time()

def checkpoint(stage):
    rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    free=shutil.disk_usage(DEST).free
    assert rss<6*1024**3, ('RSS cap',rss)
    assert free>2*1024**3, ('Free disk guard',free)
    d={'stage':stage,'elapsed_seconds':round(time.time()-started,2),'counts':dict(counts),'peak_RSS_bytes':rss,'free_disk_bytes':free}
    (OUT/'patient-extraction-progress.json').write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps(d),flush=True);stages.append(d)

assert not selected.exists() and not sorted_bam.exists(), 'Avoid overwriting an existing extraction'
with pysam.AlignmentFile(str(SRC),'rb',index_filename=str(IDX)) as bam:
    alternatives=[r for r in bam.references if ('chr6_' in r or r.startswith('HLA-'))]
    assert not alternatives, 'Unreviewed alternative contigs: '+str(alternatives)
    header=bam.header.to_dict();header['HD']['SO']='unsorted'
    with pysam.AlignmentFile(str(selected),'wb',header=header,threads=1) as dest:
        for category,iterator in [('chr6',bam.fetch('chr6')),('both_unmapped',bam.fetch('*'))]:
            for r in iterator:
                counts[category+'_input_records']+=1
                if r.flag&(256|2048):counts[category+'_secondary_or_supplementary_excluded']+=1;continue
                if category=='chr6' and (r.is_unmapped or not r.is_proper_pair):counts['chr6_not_mapped_proper_pair_excluded']+=1;continue
                if category=='both_unmapped' and (r.flag&12!=12 or not r.is_paired):counts['unmapped_not_both_unmapped_paired_excluded']+=1;continue
                if not (r.is_read1 ^ r.is_read2):counts['invalid_mate_designation_excluded']+=1;continue
                if r.query_sequence is None:counts['missing_sequence_excluded']+=1;continue
                if r.query_qualities is None:counts['missing_real_quality_excluded']+=1;continue
                assert len(r.query_sequence)==len(r.query_qualities)
                dest.write(r);counts[category+'_selected_records']+=1
                counts['selected_QC_fail_records']+=int(r.is_qcfail)
                counts['selected_duplicate_flag_records']+=int(r.is_duplicate)
                counts['selected_MAPQ_below_20_records']+=int(r.mapping_quality<20)
                if sum(counts[k] for k in ('chr6_selected_records','both_unmapped_selected_records'))%1000000==0:checkpoint('selecting')
checkpoint('selection_complete')
sort_args=['-n','-@','1','-m','512M','-T',str(DEST/'sorttmp'),'-o',str(sorted_bam),str(selected)]
pysam.sort(*sort_args)
checkpoint('name_sort_complete')
files=[DEST/'RNA_HLA.R1.fq.gz',DEST/'RNA_HLA.R2.fq.gz']
hashes=[hashlib.sha256(),hashlib.sha256()]
with pysam.AlignmentFile(str(sorted_bam),'rb') as bam, gzip.open(files[0],'wt',compresslevel=6) as f1,gzip.open(files[1],'wt',compresslevel=6) as f2:
    for name,group in itertools.groupby(bam.fetch(until_eof=True),key=lambda r:r.query_name):
        by_rg=defaultdict(lambda:defaultdict(list))
        for r in group:by_rg[r.get_tag('RG') if r.has_tag('RG') else ''][1 if r.is_read1 else 2].append(r)
        for rg,mates in by_rg.items():
            counts['RG_query_name_groups']+=1
            if any(len(x)>1 for x in mates.values()):counts['duplicate_primary_mate_groups_quarantined']+=1;continue
            if set(mates)!={1,2}:counts['orphan_RG_query_names_not_written']+=1;counts['orphan_records_not_written']+=sum(map(len,mates.values()));continue
            r1,r2=mates[1][0],mates[2][0]
            category='both_unmapped' if r1.is_unmapped and r2.is_unmapped else 'chr6'
            assert (r1.is_unmapped and r2.is_unmapped) or (not r1.is_unmapped and not r2.is_unmapped)
            name_out=rg+'|'+name if rg else name
            assert not any(ch.isspace() for ch in name_out)
            for j,(r,f) in enumerate(((r1,f1),(r2,f2))):
                seq=r.get_forward_sequence();qual=r.get_forward_qualities()
                assert qual is not None and len(seq)==len(qual)
                quality=''.join(chr(q+33) for q in qual)
                s='@'+name_out+'\n'+seq+'\n+\n'+quality+'\n'
                f.write(s);hashes[j].update(s.encode())
                counts['written_bases_R'+str(j+1)]+=len(seq)
            counts[category+'_complete_pairs_written']+=1;counts['complete_pairs_written']+=1
            if counts['complete_pairs_written']%500000==0:checkpoint('writing_paired_fastq')
checkpoint('paired_fastq_complete')
roundtrip=[hashlib.sha256(),hashlib.sha256()];n=0
with pysam.FastxFile(str(files[0])) as a, pysam.FastxFile(str(files[1])) as b:
    for r1,r2 in itertools.zip_longest(a,b):
        assert r1 is not None and r2 is not None and r1.name==r2.name
        for j,r in enumerate((r1,r2)):
            assert r.quality is not None and len(r.sequence)==len(r.quality)
            roundtrip[j].update(('@'+r.name+'\n'+r.sequence+'\n+\n'+r.quality+'\n').encode())
        n+=1
assert n==counts['complete_pairs_written']
assert [h.hexdigest() for h in hashes]==[h.hexdigest() for h in roundtrip]
assert counts['duplicate_primary_mate_groups_quarantined']==0
checkpoint('full_fastq_roundtrip_passed')
def sha256(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
result={'status':'passed','source':str(SRC),'source_index':str(IDX),'pysam_version':pysam.__version__,'samtools_version':pysam.__samtools_version__,'method':'chr6 mapped proper-paired plus both-unmapped paired records; exclude secondary/supplementary; no MAPQ threshold or duplicate/QC-fail exclusion. No missing patient qualities fabricated. Name-sort then pair by RG+query-name, quarantine duplicate mate keys; omit/count orphans. Reverse aligned sequence and qualities back to original read orientation. Full paired FASTQ name, sequence, quality roundtrip and count validated.','differences_from_official_extraction':'Primary filtering and explicit RG+name pairing/orphan auditing added. Public official extraction uses samtools, bedtools, pigz; this implementation uses pysam bundled samtools, get_forward_sequence/qualities, gzip.','reference_alternative_contigs':alternatives,'counts':dict(counts),'sort_args':sort_args,'stages':stages,'paired_FASTQ_content_sha256':[h.hexdigest() for h in roundtrip],'files':[{'path':str(p),'bytes':p.stat().st_size,'sha256':sha256(p)} for p in files],'derived_bam_paths':[str(selected),str(sorted_bam)],'derived_bam_bytes':[p.stat().st_size for p in (selected,sorted_bam)],'elapsed_seconds':time.time()-started,'limits':['Read pairs are not independently tagged molecules; no UMI deduplication.','Bulk FFPE tumor contains nonmalignant cells. HLA typing is exploratory, not clinical typing or allele retention/LOH.','Pairs with one unmapped mate or improper alignment are not in this arcas-style subset.']}
(OUT/'patient-extraction.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ('status','counts','elapsed_seconds','files')}),flush=True)
