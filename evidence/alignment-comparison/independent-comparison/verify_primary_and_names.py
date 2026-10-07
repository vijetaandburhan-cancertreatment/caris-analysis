"""Separate all-site comparison, original-name recount, and same-record PABPC1 linkage check."""
from pathlib import Path
import csv,collections,json,hashlib,datetime,time
import pysam
from audit_new_alignment import add_record,summarize_group,BITS,REVERSE,FLAGS

O=Path(__file__).resolve().parent;R=O.parent
def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def key(r):return r['chrom'],int(r['pos1'])
def js(p,x):p.write_text(json.dumps(x,indent=2))
def main():
    started=time.time();ours={key(r):r for r in rows(O/'independent-all4015-counts.tsv')};primary={key(r):r for r in rows(R/'comparison/fixed4015-comparison.tsv')}
    assert set(ours)==set(primary) and len(ours)==4015
    fieldmap={**{'new_'+c:'new_'+c for c in 'ACGT'},'new_ACGT_depth':'new_depth','new_conflict':'new_conflicting_names','new_ref':'new_ref','new_alt':'new_alt','new_other':'new_other','new_callable20':'new_callable20','new_zero_ACGT':'new_zero_ACGT'}
    diffs=[];checks=0
    for k,a in ours.items():
        b=primary[k]
        for x,y in fieldmap.items():
            checks+=1
            if int(a[x])!=int(b[y]):diffs.append({'site':k,'independent_field':x,'independent':a[x],'primary':b[y]})
        assert a['original_DNA_category']==b['DNA_category']
        if a['DNA_dominant_allele']:
            for threshold in (90,98):
                expected=int(int(a['new_callable20']) and int(a[f'new_retains{threshold}']))
                checks+=1
                if expected!=int(b[f'new_retain{threshold}_with_coverage20']):diffs.append({'site':k,'field':f'coverage20_retain{threshold}','independent':expected,'primary':b[f'new_retain{threshold}_with_coverage20']})
    assert not diffs,diffs[:10]
    dom=[r for r in ours.values() if r['DNA_dominant_allele']]
    fixed={'DNA_dominant_denominator':len(dom),'coverage_qualified_retain90':sum(int(r['new_callable20']) and int(r['new_retains90']) for r in dom),'coverage_qualified_retain98':sum(int(r['new_callable20']) and int(r['new_retains98']) for r in dom),'below20_remain_in_denominator':sum(not int(r['new_callable20']) for r in dom),'fraction_only_retain90_including_below20_secondary':sum(int(r['new_retains90']) for r in dom),'fraction_only_retain98_including_below20_secondary':sum(int(r['new_retains98']) for r in dom)}
    assert fixed=={'DNA_dominant_denominator':2873,'coverage_qualified_retain90':2805,'coverage_qualified_retain98':2721,'below20_remain_in_denominator':59,'fraction_only_retain90_including_below20_secondary':2860,'fraction_only_retain98_including_below20_secondary':2776}
    frozen=rows(R/'frozen-exception-sites.tsv');base={key(r):r for r in rows(R/'baseline-fixed-4015.tsv')};p=collections.defaultdict(dict)
    for i,r in enumerate(frozen):p[r['chrom']][int(r['pos1'])-1]=i
    oldgroups=[{} for _ in frozen];oldreadcounts=[collections.Counter() for _ in frozen]
    oldbam=R/'original-complete-names/original-3426-names.bam';newbam=R/'patient-full-pass/cohort-and-exception-records.bam'
    newqualified=json.loads((O/'independent-exception-QNAME-evidence.json').read_text());wanted=set((R/'frozen-exception-qnames.txt').read_text().splitlines());oldseen=set();oldrecords=0
    pab=('chr8',100705590,100705603);phase={}
    newcontexts=collections.defaultdict(list)
    # Independent direct aligned-pair walker; both PABPC1 substitutions must be
    # carried by this same qualifying BAM record, not merely a shared QNAME.
    for kind,bam in [('original',oldbam),('new',newbam)]:
        joint_names=set();joint_records=0;examples=[]
        with pysam.AlignmentFile(str(bam),'rb') as inp:
            refs=inp.references
            for rd in inp.fetch(until_eof=True):
                if kind=='original':
                    oldseen.add(rd.query_name);oldrecords+=1
                    if not rd.is_unmapped and refs[rd.reference_id] in p:add_record(rd,p[refs[rd.reference_id]],oldgroups,oldreadcounts)
                named=rd.query_name in wanted
                if rd.is_unmapped:
                    if kind=='new' and named:newcontexts[rd.query_name].append({'unmapped':True,'chrom':None,'start':-1,'end':-1,'base_sites':set()})
                    continue
                ch=refs[rd.reference_id];need_context=kind=='new' and named
                need_pab=ch==pab[0] and rd.reference_start<=pab[1] and rd.reference_end>pab[2]
                if not need_context and not need_pab:continue
                aligned={rp:qp for qp,rp in rd.get_aligned_pairs(matches_only=True)}
                if need_context:
                    covered={(ch,pos+1) for pos in p.get(ch,{}) if pos in aligned}
                    newcontexts[rd.query_name].append({'unmapped':False,'chrom':ch,'start':rd.reference_start,'end':rd.reference_end,'base_sites':covered})
                if not need_pab or rd.flag&FLAGS or not rd.is_proper_pair or rd.mapping_quality<30 or rd.query_sequence is None or rd.query_qualities is None:continue
                if pab[1] not in aligned or pab[2] not in aligned:continue
                q1,q2=aligned[pab[1]],aligned[pab[2]]
                if min(rd.query_qualities[q1],rd.query_qualities[q2])<25:continue
                if rd.query_sequence[q1].upper()=='G' and rd.query_sequence[q2].upper()=='C':
                    joint_names.add(rd.query_name);joint_records+=1
                    examples.append({'query_name':rd.query_name,'flag':rd.flag,'start0':rd.reference_start,'cigar':rd.cigarstring,'mapq':rd.mapping_quality,'pos591_query_index0':q1,'pos604_query_index0':q2,'BQ591':int(rd.query_qualities[q1]),'BQ604':int(rd.query_qualities[q2])})
        phase[kind]={'same_record_supporting_QNAMEs':sorted(joint_names),'same_record_supporting_records':joint_records,'record_evidence':examples}
    assert oldseen==wanted
    expectedtrans=[];statecases=collections.Counter();oldstates={}
    for i,r in enumerate(frozen):
        k=key(r);b=base[k];obs=summarize_group(oldgroups[i]);ref=b['ref'];alt=b['alt'];expected=ref if b['DNA_category']=='ref_dominant' else alt if b['DNA_category']=='alt_dominant' else None
        assert obs[ref]==int(b['RNA_fragment_ref']) and obs[alt]==int(b['RNA_fragment_alt']) and sum(obs[x] for x in 'ACGT' if x not in (ref,alt))==int(b['RNA_fragment_other'])
        states={n:REVERSE[m] if m in REVERSE else 'nonACGT' if m==0 else 'conflict' for n,m in oldgroups[i].items()};oldstates[k]=states
        newstates=newqualified[f'{k[0]}:{k[1]}']['qualifying_QNAME_states']
        for name,oldbase in states.items():
            if oldbase not in BITS:continue
            if r['review_class']=='MAPQ60_only_not_baseline_flag':
                if oldbase!=r['unexpected_base']:continue
            elif oldbase==expected:continue
            newbase=newstates.get(name)
            if newbase==oldbase:status='retains_same_unexpected_base'
            elif newbase in BITS and newbase==expected:status=f'now_{newbase}_DNA_expected'
            elif newbase in BITS:status=f'now_{newbase}_other_base'
            elif newbase=='conflict':status='conflicting_qualified_bases'
            else:
                context=newcontexts.get(name,[]);assert context
                if all(z['unmapped'] for z in context):status='only_unmapped_emitted_records'
                elif any(k in z['base_sites'] for z in context):status='site_base_overlap_but_fails_matched_filters'
                elif any(not z['unmapped'] and z['chrom']==k[0] and z['start']<=k[1]-1<z['end'] for z in context):status='site_reference_span_but_no_query_base'
                else:status='mapped_without_base_at_original_site'
            expectedtrans.append((f'{k[0]}:{k[1]}',oldbase,name,status));statecases[status]+=1
    prows=rows(R/'name-comparison/unexpected-name-transitions.tsv')
    provided=[(r['site'],r['original_unexpected_base'],r['QNAME'],r['status']) for r in prows]
    if collections.Counter(expectedtrans)!=collections.Counter(provided):
        js(O/'transition-audit-debug-differences.json',{'independent_only':list((collections.Counter(expectedtrans)-collections.Counter(provided)).elements()),'primary_only':list((collections.Counter(provided)-collections.Counter(expectedtrans)).elements())})
        raise AssertionError('Per-name transition mismatch; see independent debug differences')
    old_shared={n for n,v in oldstates[('chr8',100705591)].items() if v=='G'}&{n for n,v in oldstates[('chr8',100705604)].items() if v=='C'}
    new_shared={n for n,v in newqualified['chr8:100705591']['qualifying_QNAME_states'].items() if v=='G'}&{n for n,v in newqualified['chr8:100705604']['qualifying_QNAME_states'].items() if v=='C'}
    # Restrict same-record support to unambiguous name-level allele states.
    for kind,allowed in [('original',old_shared),('new',new_shared)]:
        phase[kind]['unambiguous_shared_QNAMEs']=sorted(allowed)
        phase[kind]['unambiguous_shared_names_with_same_record']=sorted(allowed&set(phase[kind]['same_record_supporting_QNAMEs']))
        phase[kind]['unambiguous_shared_names_without_same_record']=sorted(allowed-set(phase[kind]['same_record_supporting_QNAMEs']))
    linkage={'original_shared_names':len(old_shared),'new_shared_names':len(new_shared),'original_shared_still_both':len(old_shared&new_shared),'new_shared_not_in_old_shared':len(new_shared-old_shared),'original_shared_with_same_qualifying_record':len(phase['original']['unambiguous_shared_names_with_same_record']),'new_shared_with_same_qualifying_record':len(phase['new']['unambiguous_shared_names_with_same_record']),'new_same_record_supporting_records':phase['new']['same_record_supporting_records'],'scope':'Two linked RNA bases on the same qualifying aligned record; this does not establish tumor-specific origin, an RNA-editing mechanism, or a clinical target.'}
    ps=json.loads((R/'name-comparison/summary.json').read_text())['PABPC1_linkage'];assert linkage['original_shared_names']==ps['original_shared'] and linkage['new_shared_names']==ps['new_shared_all'] and linkage['original_shared_still_both']==ps['old_shared_still_both'] and linkage['new_shared_not_in_old_shared']==ps['new_shared_not_in_old_shared']
    js(O/'independent-PABPC1-same-record-evidence.json',phase)
    report={'status':'PASS','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'fixed_cohort_loci_compared':4015,'integer_field_comparisons':checks,'count_differences':diffs,'fixed_denominator_summary':fixed,'original_named_BAM_records':oldrecords,'original_requested_names_complete':len(oldseen),'all16_original_ACGT_vectors_match_frozen':True,'individual_old_unexpected_supporter_transitions_compared':len(expectedtrans),'all_individual_transition_classes_match_primary':True,'transition_classes':dict(statecases),'PABPC1_independent_linkage':linkage,'elapsed_seconds':round(time.time()-started,3),'source_sha256':{str(p.relative_to(R)):sha(p) for p in [Path(__file__),O/'audit_new_alignment.py',R/'comparison/fixed4015-comparison.tsv',R/'name-comparison/unexpected-name-transitions.tsv',R/'name-comparison/summary.json',oldbam,newbam]},'limits':['Only same-record linkage was independently tested; root auxiliary forward-offset/base-position tracing is a separate analysis and is not validated by this check.','Preserved baseline denominators include59 DNA-dominant sites failing new depth20. Raw fraction-only rates at shallow sites are secondary and not substituted for coverage-qualified retention.','Mapping/version/reference/parameter changes act together; no sole-D8 cause or corrected biological truth is implied.']}
    js(O/'primary-and-transition-audit.json',report);print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},indent=2))
if __name__=='__main__':main()
