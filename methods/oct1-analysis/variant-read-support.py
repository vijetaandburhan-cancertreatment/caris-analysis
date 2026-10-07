"""Targeted research read evidence; no claims of clinical validation or somatic status."""
import csv, json, pathlib, collections, statistics, time
import pysam
ROOT=pathlib.Path(__file__).resolve().parents[2]
OUT=ROOT/'work/oct1-analysis'; SOURCE=ROOT/'outputs/caris-raw-data/TN26-279853'
CAND=list(csv.DictReader((OUT/'variants.candidates.tsv').open(),delimiter='\t'))
WINDOWS={x['gene']:x for x in json.loads((OUT/'variant-read-support.hg38-reference-windows.json').read_text())}

def layout(read):
    """Reference-base -> query position; deletions and skipped intervals separate."""
    pos=read.reference_start; q=0; bases={}; dels=[]; skips=[]; insertions=[]
    for op,n in read.cigartuples or []:
        if op in (0,7,8):
            bases.update((pos+i,q+i) for i in range(n)); pos+=n; q+=n
        elif op==1: insertions.append((pos,q,n)); q+=n
        elif op==2: dels.append((pos,n)); pos+=n
        elif op==3: skips.append((pos,n)); pos+=n
        elif op==4: q+=n
        elif op in (5,6): pass
        else: raise RuntimeError(op)
    return bases,dels,skips,insertions

def setup(v):
    p=int(v['position_1based'])-1; ref=v['ref']; alt=v['alt']
    d={'p':p,'chrom':v['hg38_chrom'],'ref':ref,'alt':alt}
    if len(ref)==len(alt)==1:
        d.update({'type':'SNV','start':p,'end':p+1})
    elif len(ref)>len(alt) and ref.startswith(alt):
        win=WINDOWS[v['gene']]; seq=win['dna'].upper(); ws=win['start']; offset=p-ws
        assert seq[offset:offset+len(ref)]==ref
        deletion_start=p+len(alt); length=len(ref)-len(alt)
        expected=seq[:offset]+alt+seq[offset+len(ref):]
        eq=[s for s in range(max(ws+3,deletion_start-40),min(win['end']-length-3,deletion_start+40)+1) if seq[:s-ws]+seq[s-ws+length:]==expected]
        assert deletion_start in eq
        d.update({'type':'deletion','length':length,'equivalent_deletion_starts_0based':eq,'refwindow_start':ws,'refseq':seq,'start':min(eq)-3,'end':max(eq)+length+3})
    else: raise RuntimeError('Unsupported candidate allele '+str(v))
    return d

def classify(read,d,minq=20):
    if read.is_unmapped or read.query_sequence is None: return 'unassessable',None
    bases,dels,skips,insertions=layout(read); seq=read.query_sequence; qual=read.query_qualities
    def qok(positions): return qual is not None and all(qual[q]>=minq for q in positions)
    if d['type']=='SNV':
        p=d['p']; q=bases.get(p)
        if q is None:
            return ('deletion_at_locus' if any(s<=p<s+n for s,n in dels) else 'splice_skip' if any(s<=p<s+n for s,n in skips) else 'unassessable'),None
        if not qok([q]): return 'low_base_quality',q
        allele='alt' if seq[q]==d['alt'] else 'ref' if seq[q]==d['ref'] else 'other_base'
        return allele,q
    # Deletion candidates: use all locally equivalent representations, with three
    # reference-matching, Q>=20 aligned bases on either side of the deletion.
    for s,n in dels:
        if n==d['length'] and s in d['equivalent_deletion_starts_0based']:
            flank=list(range(s-3,s))+list(range(s+n,s+n+3))
            if not all(p in bases for p in flank): return 'indel_insufficient_flank',None
            qs=[bases[p] for p in flank]
            if not qok(qs): return 'low_base_quality',None
            if any(seq[bases[p]] != d['refseq'][p-d['refwindow_start']] for p in flank): return 'indel_flank_mismatch',None
            return 'alt',bases[s-1]
    region=list(range(d['start'],d['end']))
    if any(s<d['end'] and s+n>d['start'] for s,n in skips): return 'splice_skip',None
    if any(s<d['end'] and s+n>d['start'] for s,n in dels): return 'other_indel',None
    if any(d['start']<s<d['end'] for s,q,n in insertions): return 'other_indel',None
    if not all(p in bases for p in region): return 'indel_insufficient_flank',None
    qs=[bases[p] for p in region]
    if not qok(qs): return 'low_base_quality',None
    if any(seq[bases[p]]!=d['refseq'][p-d['refwindow_start']] for p in region): return 'indel_region_mismatch',None
    return 'ref',bases[d['p']]

def summarize(bam,v,d,kind):
    raw=collections.Counter(); exclusions=collections.Counter(); eligible=collections.Counter()
    orient=collections.defaultdict(collections.Counter); frag=collections.defaultdict(list)
    detail=collections.defaultdict(list); clean_frag=collections.defaultdict(set); total=0
    for r in bam.fetch(d['chrom'],d['start'],d['end']):
        total+=1
        raw_call,_=classify(r,d,0); raw[raw_call]+=1
        reason=('unmapped' if r.is_unmapped else 'secondary' if r.is_secondary else 'supplementary' if r.is_supplementary else 'QC_failed' if r.is_qcfail else 'duplicate' if r.is_duplicate else 'MAPQ_below_20' if r.mapping_quality<20 else None)
        if reason: exclusions[reason]+=1; continue
        call,q=classify(r,d,20); eligible[call]+=1
        if call in ('alt','ref','other_base','other_indel','indel_region_mismatch','indel_flank_mismatch'):
            rg=r.get_tag('RG') if r.has_tag('RG') else ''
            frag[(rg,r.query_name)].append(call)
        if call in ('alt','ref'):
            orient[call]['reverse' if r.is_reverse else 'forward']+=1
            soft=any(op==4 for op,n in r.cigartuples or [])
            near= q is not None and min(q-r.query_alignment_start,r.query_alignment_end-1-q)<5
            if not soft and not near:
                clean_frag[call].add((r.get_tag('RG') if r.has_tag('RG') else '',r.query_name))
            detail[call].append({'softclipped':soft,'near_aligned_end_lt5bp':near,'mapq':r.mapping_quality,'baseq':int(r.query_qualities[q]) if q is not None else None})
    fc=collections.Counter(); discordant=0
    for calls in frag.values():
        choices=set(calls)
        if len(choices)==1: fc[calls[0]]+=1
        else: fc['discordant']+=1; discordant+=1
    ed=eligible['alt']+eligible['ref']; fd=fc['alt']+fc['ref']
    extras={}
    for allele in ('ref','alt'):
        rows=detail[allele]
        extras[allele]={'forward_reads':orient[allele]['forward'],'reverse_reads':orient[allele]['reverse'],'softclipped_reads':sum(x['softclipped'] for x in rows),'support_near_read_end_lt5bp':sum(x['near_aligned_end_lt5bp'] for x in rows),'median_mapq':statistics.median(x['mapq'] for x in rows) if rows else None,'median_support_baseq':statistics.median(x['baseq'] for x in rows if x['baseq'] is not None) if rows else None,'fragments_with_a_supporting_read_without_softclip_or_near_end':len(clean_frag[allele])}
    return {'kind':kind,'gene':v['gene'],'protein':v['protein'],'chrom':d['chrom'],'position_1based':d['p']+1,'ref':d['ref'],'alt':d['alt'],'variant_type':d['type'],'fetch_start_0based':d['start'],'fetch_end_0based_exclusive':d['end'],'equivalent_deletion_starts_0based':d.get('equivalent_deletion_starts_0based'),'raw_alignments_fetched':total,'raw_allele_counts_no_flag_mapq_baseq_filter':dict(raw),'read_exclusions_first_matching_reason':dict(exclusions),'filtered_read_counts':dict(eligible),'filtered_fragment_counts':dict(fc),'assessable_ref_alt_reads':ed,'alt_read_fraction_ref_alt_only':eligible['alt']/ed if ed else None,'assessable_ref_alt_fragments':fd,'alt_fragment_fraction_ref_alt_only':fc['alt']/fd if fd else None,'allele_detail':extras,'caris_reported_DNA_vaf':float(v['tumor_vaf']),'caris_gene_TPM':float(v['gene_TPM'])}

def synthetic_checks():
    # Meaningful coordinate/CIGAR guard: BAP1/RASA1 references, exact and shifted
    # deletion representations, SNV and low-base quality must classify as expected.
    for v in CAND:
        if v['gene'] not in ('BAP1','RASA1','APC'): continue
        d=setup(v)
        if d['type']=='SNV':
            for base,result in ((d['ref'],'ref'),(d['alt'],'alt')):
                r=pysam.AlignedSegment();r.query_name='test';r.reference_start=d['p'];r.cigarstring='1M';r.query_sequence=base;r.query_qualities=[30]
                assert classify(r,d)==(result,0)
                r.query_qualities=[10]; assert classify(r,d)[0]=='low_base_quality'
        else:
            for s in d['equivalent_deletion_starts_0based']:
                st=d['start']-5;en=d['end']+5; ws=d['refwindow_start'];refseq=d['refseq'][st-ws:en-ws]
                r=pysam.AlignedSegment();r.query_name='test';r.reference_start=st;r.cigarstring=f'{s-st}M{d["length"]}D{en-s-d["length"]}M';r.query_sequence=refseq[:s-st]+refseq[s-st+d['length']:];r.query_qualities=[30]*len(r.query_sequence)
                assert classify(r,d)[0]=='alt',(v['gene'],s,classify(r,d))
            r=pysam.AlignedSegment();r.query_name='test';r.reference_start=st;r.cigarstring=f'{en-st}M';r.query_sequence=refseq;r.query_qualities=[30]*len(refseq)
            assert classify(r,d)[0]=='ref'

synthetic_checks()
results={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'method':{'reference_build':'hg38','reference_windows_source':'UCSC sequence API, two bounded public reference windows only','min_mapq':20,'min_baseq':20,'excluded_flags':['unmapped','secondary','supplementary','QC failed','duplicate'],'fragments':'group qualifying evidence by read group + query name; discordant mates excluded from ref+alt fraction','indels':'CIGAR deletion with locally equivalent shift positions and 3bp high-quality matching flanks; reference requires continuous high-quality matching interval across union of equivalent positions','RNA_warning':'RNA BAM has no duplicate flags; fragment deduplication merges mates, not PCR duplicates/UMIs. RNA library strand can bias read orientation. Low coverage/no alternate support is not proof of absent expression.','fraction_warning':'Research strict-filter fractions can differ from Caris due to filters, local indel handling and mate collapsing. Denominators exclude nonref/nonalt/ambiguous cases.','synthetic_coordinate_and_indel_checks':'passed'},'results':[]}
for kind in ('RNA','DNA'):
    path=SOURCE/f'{kind}_TN26-279853.bam'; idx=OUT/f'{path.name}.bai'
    with pysam.AlignmentFile(str(path),'rb',index_filename=str(idx),threads=2) as bam:
        for v in CAND:
            result=summarize(bam,v,setup(v),kind);results['results'].append(result)
            print(kind,v['gene'],v['protein'],result['filtered_read_counts'],result['filtered_fragment_counts'],flush=True)
(OUT/'variant-read-support.json').write_text(json.dumps(results,indent=2)+'\n')
fields=['kind','gene','protein','chrom','position_1based','ref','alt','raw_alignments_fetched','ref_reads','alt_reads','assessable_ref_alt_reads','alt_read_fraction_ref_alt_only','ref_fragments','alt_fragments','discordant_fragments','assessable_ref_alt_fragments','alt_fragment_fraction_ref_alt_only','alt_forward_reads','alt_reverse_reads','alt_softclipped_reads','alt_near_aligned_end_lt5bp','caris_reported_DNA_vaf','caris_gene_TPM']
with (OUT/'variant-read-support.tsv').open('w') as f:
    writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t');writer.writeheader()
    for r in results['results']:
        row={k:r.get(k) for k in fields};counts=r['filtered_read_counts'];frags=r['filtered_fragment_counts'];alt=r['allele_detail']['alt']
        row.update(ref_reads=counts.get('ref',0),alt_reads=counts.get('alt',0),ref_fragments=frags.get('ref',0),alt_fragments=frags.get('alt',0),discordant_fragments=frags.get('discordant',0),alt_forward_reads=alt['forward_reads'],alt_reverse_reads=alt['reverse_reads'],alt_softclipped_reads=alt['softclipped_reads'],alt_near_aligned_end_lt5bp=alt['support_near_read_end_lt5bp']);writer.writerow(row)
