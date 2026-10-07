"""Independent bounded recount using pysam aligned-pairs/haplotype extraction.
Does not import or call the audited analysis script. Reads three loci in each BAM.
"""
import json, pathlib, collections, hashlib
import pysam

out=pathlib.Path(__file__).resolve().parent
source=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
reported=json.loads((out/'variant-read-support.json').read_text())
candidates={r['gene']:r for r in json.loads((out/'variants.candidates.json').read_text()) if r['gene'] in ('APC','BAP1','RASA1')}
windows={r['gene']:r for r in json.loads((out/'variant-read-support.hg38-reference-windows.json').read_text())}

def allele_read(read, v):
    pos=v['position_1based']-1
    # Pysam query_sequence is stored in reference-aligned orientation in BAM,
    # including reverse-strand alignments; no second reverse-complement.
    pairs=read.get_aligned_pairs()
    refmap={r:q for q,r in pairs if r is not None and q is not None}
    if len(v['ref'])==len(v['alt'])==1:
        if pos not in refmap:return 'unassessable',None
        q=refmap[pos]
        if read.query_qualities[q]<20:return 'low_base_quality',None
        base=read.query_sequence[q]
        return ('alt' if base==v['alt'] else 'ref' if base==v['ref'] else 'other_base'), None
    # Haplotype across the same three reference flanks. A retained local
    # insertion makes the sequence differ and is not counted as simple deletion.
    start=pos+len(v['alt'])
    end=pos+len(v['ref'])
    lo=start-3;hi=end+3
    if lo not in refmap or hi-1 not in refmap:return 'insufficient_flank',None
    qlo,qhi=refmap[lo],refmap[hi-1]+1
    if qhi<qlo:return 'bad_query_order',None
    refseq=windows[v['gene']]['dna'][lo-windows[v['gene']]['start']:hi-windows[v['gene']]['start']]
    assert len(refseq)==hi-lo
    expected=refseq[:3]+refseq[-3:]
    obs=read.query_sequence[qlo:qhi]
    qualities=read.query_qualities[qlo:qhi]
    # Track only actual CIGAR deletions vs splice skips, independently from
    # aligned-pairs where both produce absent query positions.
    rp=read.reference_start;del_intervals=[];skip_intervals=[]
    for op,n in read.cigartuples:
        if op==2:del_intervals.append((rp,rp+n))
        if op==3:skip_intervals.append((rp,rp+n))
        if op in (0,2,3,7,8):rp+=n
    if any(a<hi and b>lo for a,b in skip_intervals):return 'splice_skip',None
    if min(qualities)<20:return 'low_base_quality',None
    if obs==expected and (start,end) in del_intervals:return 'alt',None
    if obs==refseq and not any(a<hi and b>lo for a,b in del_intervals):return 'ref',None
    return 'other_local_haplotype', {'observed_length':len(obs),'expected_ref_length':len(refseq),'expected_alt_length':len(expected)}

records=[]
for kind in ('RNA','DNA'):
    with pysam.AlignmentFile(str(source/f'{kind}_TN26-279853.bam'),'rb',index_filename=str(out/f'{kind}_TN26-279853.bam.bai')) as b:
        for g,v in candidates.items():
            old=next(r for r in reported['results'] if r['gene']==g and r['kind']==kind)
            counts=collections.Counter();frags=collections.defaultdict(set);flagcounts=collections.Counter();nh=collections.Counter();mq=collections.Counter();orient=collections.Counter();cigar=collections.Counter();proper=collections.Counter();supportnames=collections.defaultdict(set)
            readcoords=collections.defaultdict(collections.Counter);fragcoords=collections.defaultdict(lambda:collections.defaultdict(set));strictfrags=collections.defaultdict(set);strictfragcoords=collections.defaultdict(lambda:collections.defaultdict(set))
            for read in b.fetch(v['hg38_chrom'],old['fetch_start_0based'],old['fetch_end_0based_exclusive']):
                if read.flag & (4|256|512|1024|2048):
                    flagcounts['excluded_flags']+=1;continue
                if read.mapping_quality<20:
                    flagcounts['low_mapping_quality']+=1;continue
                call,detail=allele_read(read,v);counts[call]+=1
                if call in ('ref','alt','other_base','other_local_haplotype'):
                    key=(read.get_tag('RG') if read.has_tag('RG') else '',read.query_name)
                    frags[key].add(call)
                if call in ('ref','alt'):
                    nh[str(read.get_tag('NH')) if read.has_tag('NH') else 'missing']+=1
                    mq[read.mapping_quality]+=1
                    orient[(call,'rev' if read.is_reverse else 'fwd')]+=1
                    proper[(call,read.is_proper_pair)]+=1
                    if call=='alt':cigar[read.cigarstring]+=1
                    supportnames[call].add(read.query_name)
                    readcoords[call][(read.reference_start,read.reference_end,read.cigarstring,read.is_reverse,read.is_read1)]+=1
                    coord=(read.reference_id,min(read.reference_start,read.next_reference_start),abs(read.template_length),read.is_reverse if read.is_read1 else not read.is_reverse) if read.is_paired and read.next_reference_id==read.reference_id and not read.mate_is_unmapped else (read.reference_id,read.reference_start,read.reference_end,read.is_reverse)
                    fragcoords[call][coord].add(read.query_name)
                    rm={rp:qp for qp,rp in read.get_aligned_pairs() if qp is not None and rp is not None}
                    qp=rm.get(v['position_1based']-1)
                    soft=any(op==4 for op,n in read.cigartuples)
                    near=qp is not None and min(qp-read.query_alignment_start,read.query_alignment_end-1-qp)<5
                    if not soft and not near:
                        strictfrags[call].add(read.query_name)
                        strictfragcoords[call][coord].add(read.query_name)
            fc=collections.Counter(next(iter(cs)) if len(cs)==1 else 'discordant' for cs in frags.values())
            comparable={a:(counts[a],old['filtered_read_counts'].get(a,0)) for a in ('ref','alt')}
            comparable_frag={a:(fc[a],old['filtered_fragment_counts'].get(a,0)) for a in ('ref','alt')}
            assert all(a==b for a,b in comparable.values()),(g,kind,comparable)
            assert all(a==b for a,b in comparable_frag.values()),(g,kind,comparable_frag)
            records.append({'kind':kind,'gene':g,'independent_read_counts':dict(counts),'independent_fragment_counts':dict(fc),'read_counts_match_original':True,'fragment_counts_match_original':True,'independent_support_mapping_quality_counts':dict(mq),'independent_support_NH_counts':dict(nh),'independent_support_read_orientation':{str(k):val for k,val in orient.items()},'support_proper_pair':{str(k):val for k,val in proper.items()},'alternate_CIGAR_examples':cigar.most_common(5),'flags':dict(flagcounts),'coordinate_diversity':{a:{'distinct_read_start_end_cigar_orientation_readnumber':len(readcoords[a]),'distinct_inferred_fragment_coordinate_groups':len(fragcoords[a]),'largest_inferred_fragment_coordinate_group_distinct_querynames':max(map(len,fragcoords[a].values()),default=0),'strict_no_softclip_or_end_support_querynames':len(strictfrags[a]),'strict_distinct_inferred_fragment_coordinate_groups':len(strictfragcoords[a])} for a in ('ref','alt')},'coordinate_diversity_definition':'Read key=start,end,CIGAR,orientation,read1 flag. Fragment surrogate=chrom,min(read start,mate start),abs(template length),read1 orientation; not UMI deduplication, independent molecule count or transcript-aware PCR family certainty. Strict flag means a supporting read without softclip and variant anchor >=5bp from aligned ends.'})

q=json.loads((out/'bam-qc.json').read_text())
flagstat_audit=[]
for f in q['files']:
    x=f['flagstat']['QC-passed reads'];idx=f['idxstats']
    flagstat_audit.append({'kind':f['kind'],'total_equals_primary_plus_secondary_plus_supplementary':x['total']==x['primary']+x['secondary']+x['supplementary'],'primary_equals_read1_plus_read2':x['primary']==x['read1']+x['read2'],'indexed_mapped_matches_flagstat':sum(a['mapped'] for a in idx)==x['mapped'],'primary_read1_count':x['read1'],'mapped_primary_percent':x['primary mapped %'],'duplicate_flagged':x['duplicates'],'warning':'RNA zero duplicate flags means unmarked duplicates, not zero PCR duplication' if f['kind']=='RNA' else 'Duplicate-flag exclusions respected in target recount'})

result={'status':'bounded independent recount passed','audited_script_sha256':hashlib.sha256((out/'variant-read-support.py').read_bytes()).hexdigest(),'independent_method':'pysam aligned-pairs mapping and full local reference/alternate haplotype string, not imported original functions; baseQ20, mapQ20, primary QC-passing nonduplicate reads; RG+query name collapse','results':records,'bam_qc_arithmetic_audit':flagstat_audit,'limitations':['No matched normal: RNA support does not prove tumor-specific or somatic status.','Mate collapse is not UMI/PCR molecule deduplication, especially RNA with no duplicate flags.','MAPQ255 is STAR unique mapping convention; not a literal phred confidence of255. NH examined independently at three loci.','Ref+alt fractions exclude reads failing filters or unable to span the required region and are not cancer-cell fractions or population clonality estimates.','BAP1 deletion fraction can have reference/alternate span ascertainment bias; counts establish expression support, not absolute RNA allele fraction.','Alignment support does not establish full-length functional transcript, translated protein, MHC presentation, immunogenicity or vaccine effectiveness.','BAM quickcheck/index/flagstat support file usability but do not establish original library quality, capture coverage, purity, absence of sample mixup or complete exome representation.','The all21 RNA-positive statement is script/result checked, but only three loci independently recounted here.']}
allrna=[r for r in reported['results'] if r['kind']=='RNA']
result['all21_RNA_positive_script_result_audit']={'record_count':len(allrna),'all_positive_alt_read_counts':all(r['filtered_read_counts'].get('alt',0)>0 for r in allrna),'all_positive_alt_fragment_counts':all(r['filtered_fragment_counts'].get('alt',0)>0 for r in allrna),'minimum_alt_fragment_count':min(r['filtered_fragment_counts'].get('alt',0) for r in allrna),'counts_not_independently_recomputed_for_other18':True}
(out/'variants.peer-review.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
