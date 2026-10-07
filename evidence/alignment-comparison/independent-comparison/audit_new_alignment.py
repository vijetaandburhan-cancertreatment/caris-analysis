"""Independent fixed-cohort RNA recount using aligned query/reference pairs, not CIGAR blocks."""
from pathlib import Path
import collections,csv,datetime,hashlib,json,resource,sys,time
import pysam

O=Path(__file__).resolve().parent;R=O.parent;C=R.parent
FLAGS=0xF0C;BITS={'A':1,'C':2,'G':4,'T':8};REVERSE={v:k for k,v in BITS.items()}
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def js(p,x):p.write_text(json.dumps(x,indent=2))
def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def tsv(p,rs,fields=None):
    fields=fields or (list(rs[0]) if rs else [])
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def add_record(rd,positions,groups,reads):
    if rd.flag&FLAGS or not rd.is_proper_pair or rd.mapping_quality<30 or rd.query_sequence is None or rd.query_qualities is None:return
    seq=rd.query_sequence;qual=rd.query_qualities;name=rd.query_name
    for qp,rp in rd.get_aligned_pairs(matches_only=True):
        if rp not in positions or qual[qp]<25:continue
        i=positions[rp];base=seq[qp].upper();groups[i][name]=groups[i].get(name,0)|BITS.get(base,0);reads[i][base]+=1
def summarize_group(g):
    out={x:0 for x in ['A','C','G','T','conflict','nonACGT']}
    for mask in g.values():out[REVERSE[mask] if mask in REVERSE else 'nonACGT' if mask==0 else 'conflict']+=1
    return out
def qualify():
    p=C/'oct5-population-panel-concordance-v1/qualification/synthetic-control/control.bam';groups=[{}];reads=[collections.Counter()]
    with pysam.AlignmentFile(str(p),'rb') as b:
        for rd in b:add_record(rd,{100:0},groups,reads)
    obs=summarize_group(groups[0]);assert obs=={'A':1,'C':0,'G':1,'T':0,'conflict':1,'nonACGT':0};assert dict(reads[0])=={'A':3,'G':3}
    q={'status':'PASS','method':'get_aligned_pairs(matches_only=True), independent bit-mask query-name collapse','control_bam_sha256':sha(p),'fragment_counts':obs,'read_base_counts':dict(reads[0]),'excludes':'Unmapped/mate-unmapped/secondary/QCfail/duplicate/supplementary, improper pair, MAPQ<30 or BQ<25','patient_reads_accessed':False}
    js(O/'counter-qualification.json',q);return q
def run():
    started=time.time();runpath=R/'patient-full-pass/run.json';run=json.loads(runpath.read_text());assert run['status']=='COMPLETE'
    qualify();base=rows(R/'baseline-fixed-4015.tsv');assert len(base)==4015
    assert len({(r['chrom'],int(r['pos1'])) for r in base})==4015
    wanted=set((R/'frozen-exception-qnames.txt').read_text().splitlines());assert len(wanted)==3426
    exceptions={(r['chrom'],int(r['pos1'])) for r in rows(R/'frozen-exception-sites.tsv')};assert len(exceptions)==16
    pos=collections.defaultdict(dict)
    for i,r in enumerate(base):pos[r['chrom']][int(r['pos1'])-1]=i
    groups=[{} for _ in base];readcounts=[collections.Counter() for _ in base];name_records=collections.defaultdict(list)
    source=R/'patient-full-pass/cohort-and-exception-records.bam';pysam.quickcheck(str(source));stream_counts=collections.Counter();headers={}
    with pysam.AlignmentFile(str(source),'rb') as bam:
        headers=bam.header.to_dict();refs=bam.references
        for rd in bam.fetch(until_eof=True):
            stream_counts['all_captured_records']+=1
            if rd.query_name in wanted:
                name_records[rd.query_name].append({'flag':rd.flag,'read1':rd.is_read1,'read2':rd.is_read2,'unmapped':rd.is_unmapped,'secondary':rd.is_secondary,'supplementary':rd.is_supplementary,'reference':refs[rd.reference_id] if rd.reference_id>=0 else '*','start0':rd.reference_start,'cigar':rd.cigarstring,'mapq':rd.mapping_quality,'mate_reference':refs[rd.next_reference_id] if rd.next_reference_id>=0 else '*','mate_start0':rd.next_reference_start,'template_length':rd.template_length,'record_SAM_sha256':hashlib.sha256(rd.to_string().encode()).hexdigest()})
            if rd.is_unmapped:continue
            ch=refs[rd.reference_id]
            if ch not in pos:continue
            add_record(rd,pos[ch],groups,readcounts)
    assert stream_counts['all_captured_records']==run['capture_summary']['counts']['captured_records']
    assert sum(len(v) for v in name_records.values())==run['capture_summary']['counts']['named_records']
    assert len(name_records)==run['capture_summary']['exception_names_seen']
    independently_named_flags=collections.Counter(str(record['flag']) for rs in name_records.values() for record in rs)
    independently_named_locations=collections.Counter('unmapped' if record['unmapped'] else record['reference'] for rs in name_records.values() for record in rs)
    assert dict(independently_named_flags)==run['capture_summary']['named_record_flags']
    assert dict(independently_named_locations)==run['capture_summary']['named_record_locations']
    newrows=[];fraction_summary=collections.Counter();details={}
    for i,r in enumerate(base):
        n=summarize_group(groups[i]);ref=r['ref'];alt=r['alt'];depth=sum(n[c] for c in 'ACGT');category=r['DNA_category']
        olddepth=int(r['RNA_ACGT_depth']);row={'chrom':r['chrom'],'pos1':int(r['pos1']),'ref':ref,'alt':alt,'original_DNA_category':category,'original_RNA_ACGT_depth':olddepth,'original_RNA_ref':int(r['RNA_fragment_ref']),'original_RNA_alt':int(r['RNA_fragment_alt']),'original_RNA_other':int(r['RNA_fragment_other'])}
        row.update({'new_'+c:n[c] for c in n});row['new_ref']=n[ref];row['new_alt']=n[alt];row['new_other']=sum(n[c] for c in 'ACGT' if c not in (ref,alt));row['new_ACGT_depth']=depth;row['new_callable20']=int(depth>=20);row['new_zero_ACGT']=int(depth==0);row['fixed_exception_site']=int((r['chrom'],int(r['pos1'])) in exceptions)
        row['same_ref_alt_other_counts']=int(row['new_ref']==row['original_RNA_ref'] and row['new_alt']==row['original_RNA_alt'] and row['new_other']==row['original_RNA_other'])
        row['DNA_dominant_allele']='';row['new_DNA_dominant_allele_fraction']='';row['new_retains90']='';row['new_retains98']=''
        if category in ('ref_dominant','alt_dominant'):
            allele=ref if category=='ref_dominant' else alt;row['DNA_dominant_allele']=allele
            if depth:
                f=n[allele]/depth;row['new_DNA_dominant_allele_fraction']=f;row['new_retains90']=int(f>=.90);row['new_retains98']=int(f>=.98)
                fraction_summary['nonzero_dominant_sites']+=1;fraction_summary['retains90_nonzero']+=int(f>=.90);fraction_summary['retains98_nonzero']+=int(f>=.98)
                if depth>=20:fraction_summary['callable20_dominant_sites']+=1;fraction_summary['retains90_callable20']+=int(f>=.90);fraction_summary['retains98_callable20']+=int(f>=.98)
            else:fraction_summary['zero_dominant_sites']+=1
        newrows.append(row)
        if row['fixed_exception_site']:
            details[f"{r['chrom']}:{r['pos1']}"]={'original_fixed_row':r,'new_counts':row,'qualifying_QNAME_states':{name:REVERSE[m] if m in REVERSE else 'nonACGT' if m==0 else 'conflict' for name,m in sorted(groups[i].items())},'new_read_base_counts':dict(readcounts[i])}
    assert len(newrows)==4015 and len(details)==16
    tsv(O/'independent-all4015-counts.tsv',newrows);tsv(O/'independent-16-exception-counts.tsv',[r for r in newrows if r['fixed_exception_site']]);js(O/'independent-exception-QNAME-evidence.json',details)
    inv=[]
    for name in sorted(wanted):
        rs=name_records.get(name,[]);r1=any(r['read1'] for r in rs);r2=any(r['read2'] for r in rs)
        inv.append({'query_name':name,'emitted_records_captured':len(rs),'read1_seen':r1,'read2_seen':r2,'both_mate_flags_seen':r1 and r2,'unmapped_records':sum(r['unmapped'] for r in rs),'secondary_records':sum(r['secondary'] for r in rs),'supplementary_records':sum(r['supplementary'] for r in rs),'flags':';'.join(str(x) for x in sorted({r['flag'] for r in rs})),'references':';'.join(sorted({r['reference'] for r in rs}))})
    tsv(O/'independent-frozen-QNAME-inventory.tsv',inv);js(O/'independent-frozen-QNAME-record-fields.json',dict(name_records))
    missing=sorted(wanted-set(name_records));primarymissing=run['capture_summary']['missing_exception_names'];assert missing==primarymissing
    result={'status':'INDEPENDENT_COUNT_COMPLETE_PENDING_PRIMARY_COMPARISON','fixed_cohort_loci':4015,'independent_output_rows':len(newrows),'original_DNA_categories_fixed':dict(collections.Counter(r['original_DNA_category'] for r in newrows)),'new_callable20_sites':sum(r['new_callable20'] for r in newrows),'new_zero_ACGT_sites':sum(r['new_zero_ACGT'] for r in newrows),'original_and_new_ref_alt_other_vectors_equal':sum(r['same_ref_alt_other_counts'] for r in newrows),'fixed_denominator_retention':dict(fraction_summary),'captured_records_read':stream_counts['all_captured_records'],'requested_exception_names':len(wanted),'observed_exception_names':len(name_records),'missing_exception_names':missing,'exception_names_with_both_mate_flags':sum(r['both_mate_flags_seen'] for r in inv),'all16_exception_rows_retained':len(details)==16,'runtime_seconds':round(time.time()-started,3),'peak_RSS_bytes_macOS':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'limits':['Counts are A/C/G/T query-name units, not UMI-independent molecules.','Fixed baseline4015 sites remain the denominator; new depth losses and zero counts are retained. Fractions among nonzero or callable subsets are separately labelled.','Ordinary D/N-only alignment records were not captured and nonACGT counters therefore do not represent complete whole-BAM skip/deletion inventories.','Observed requested names/mates describe records emitted by STAR and preserved by capture, not all possible mappings or specimen identity.','This is one RNA library processed through different mapping pipelines, not an independent sample. No correction or D8-only attribution follows.']}
    js(O/'independent-summary.json',result)
    js(O/'provenance.json',{'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha(Path(__file__)),'run_json_sha256':sha(runpath),'design_sha256':sha(R/'design-frozen.json'),'cohort_sha256':sha(R/'baseline-fixed-4015.tsv'),'frozen_qnames_sha256':sha(R/'frozen-exception-qnames.txt'),'BAM_sha256':sha(source),'captured_BAM_header':headers,'qualification_sha256':sha(O/'counter-qualification.json'),'outputs':{p.name:sha(p) for p in O.glob('*') if p.is_file() and p.suffix in ('.json','.tsv') and p.name not in ('provenance.json',)}})
    print(json.dumps(result,indent=2))
if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='qualify':print(json.dumps(qualify(),indent=2))
    else:run()
