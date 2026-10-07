"""Bounded alignment-dependence check; execute patient reads only after primary count completion."""
from pathlib import Path
import collections,csv,datetime,hashlib,json,resource,sys,time
import pysam

OUT=Path(__file__).resolve().parent;ROOT=OUT.parent
CACHE=Path('/Users/burhanazeem/.local/share/codex/caris-analysis')
WORK=Path('/Users/burhanazeem/Documents/Codex/2026-09-05/finances-plugin-finances-openai-curated-remote-3')
OLD=CACHE/'TN26-279853/RNA_TN26-279853.bam'
OLD_INDEX=WORK/'work/oct1-analysis/RNA_TN26-279853.bam.bai'
NEW=CACHE/'oct4-followup/genome-fusion/patient-D8/target-genes.sorted.bam'
REGIONS=[('chr3',52401007,52410008,'BAP1'),('chr5',87267882,87391931,'RASA1')]
EXCLUDED=0xF0C
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write_json(path,obj):path.write_text(json.dumps(obj,indent=2))
def within(row):
    return [gene for chrom,start0,end0,gene in REGIONS if row['chrom']==chrom and start0<=int(row['pos1'])-1<end0]
def qpos_at(read,pos0):
    r=read.reference_start;q=0
    for op,n in read.cigartuples or []:
        if op in (0,7,8):
            if r<=pos0<r+n:return q+(pos0-r),'aligned_base'
            r+=n;q+=n
        elif op in (2,3):
            if r<=pos0<r+n:return None,'deletion' if op==2 else 'reference_skip'
            r+=n
        elif op in (1,4):q+=n
        elif op in (5,6):pass
        else:raise ValueError('Unexpected CIGAR operation '+str(op))
    return None,'no_query_base'
def count_one(bam,chrom,pos1):
    all_by_name=collections.defaultdict(list);eligible=collections.defaultdict(list);reasons=collections.Counter()
    for r in bam.fetch(chrom,pos1-1,pos1):
        qp,location=qpos_at(r,pos1-1)
        base=r.query_sequence[qp].upper() if qp is not None and r.query_sequence else None
        bq=int(r.query_qualities[qp]) if qp is not None and r.query_qualities is not None else None
        exclusions=[]
        if r.flag&EXCLUDED:exclusions.append('excluded_flag')
        if not r.is_proper_pair:exclusions.append('not_proper_pair')
        if r.mapping_quality<30:exclusions.append('MAPQ_below30')
        if qp is None:exclusions.append(location)
        elif bq is None or bq<25:exclusions.append('BQ_below25_or_absent')
        d={'flag':r.flag,'start0':r.reference_start,'cigar':r.cigarstring,'mapq':r.mapping_quality,'query_position0':qp,'base':base,'BQ':bq,'excluded_reasons':exclusions}
        all_by_name[r.query_name].append(d)
        reasons.update(exclusions)
        if not exclusions:eligible[r.query_name].append(base)
    names={};counts=collections.Counter({x:0 for x in ['A','C','G','T','conflict','nonACGT']});read_counts=collections.Counter()
    for name,bases in eligible.items():
        read_counts.update(bases);a={b for b in bases if b in 'ACGT'}
        state='conflict' if len(a)>1 else next(iter(a)) if len(a)==1 else 'nonACGT'
        names[name]=state;counts[state]+=1
    return {'fragment_counts':dict(counts),'read_base_counts':dict(read_counts),'qname_states':names,'overlap_record_reasons':dict(reasons),'all_overlapping_records_by_name':dict(all_by_name)}
def qualify():
    p=ROOT/'qualification/synthetic-control/control.bam'
    with pysam.AlignmentFile(str(p),'rb') as b:r=count_one(b,'chrTest',101)
    assert r['fragment_counts']=={'A':1,'C':0,'G':1,'T':0,'conflict':1,'nonACGT':0},r['fragment_counts']
    assert r['read_base_counts']=={'A':3,'G':3}
    q={'status':'PASS','control_bam_sha256':sha(p),'fragment_counts':r['fragment_counts'],'read_base_counts':r['read_base_counts'],'qualification':'Independent CIGAR walker matches synthetic mate-collapse, conflict and flag/MAPQ/BQ expectations; no patient read access.'}
    write_json(OUT/'qualification.json',q);return q
def run():
    started=time.time();complete=json.loads((ROOT/'count-complete.json').read_text());assert complete['status']=='PASS'
    assert (ROOT/'joint-callable.tsv').exists();qualify()
    with (ROOT/'input/public-common-SNP-panel.tsv').open() as f:universe=[r for r in csv.DictReader(f,delimiter='\t') if within(r)]
    with (ROOT/'joint-callable.tsv').open() as f:selected=[r for r in csv.DictReader(f,delimiter='\t') if within(r)]
    assert len(universe)==73
    coverage_fields=['chrom','pos1','ref','alt','scope','DNA_ACGT_depth','RNA_ACGT_depth','DNA_callable20','RNA_callable20','joint_callable']
    with (ROOT/'all-public-loci-counts.tsv').open() as f:
        target_coverage=[{**{k:r[k] for k in coverage_fields},'target_gene':';'.join(within(r))} for r in csv.DictReader(f,delimiter='\t') if within(r)]
    assert len(target_coverage)==73
    with (OUT/'target-panel-coverage.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=coverage_fields+['target_gene'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(target_coverage)
    selected.sort(key=lambda r:(r['chrom'],int(r['pos1'])))
    results=[];details={};reference_checks={}
    with pysam.AlignmentFile(str(OLD),'rb',index_filename=str(OLD_INDEX)) as a,pysam.AlignmentFile(str(NEW),'rb') as b:
        for chrom,_,_,_ in REGIONS:
            reference_checks[chrom]={'original_length':a.get_reference_length(chrom),'new_length':b.get_reference_length(chrom)}
            assert a.get_reference_length(chrom)==b.get_reference_length(chrom)
        headers={'original':a.header.to_dict(),'new':b.header.to_dict()}
        for site in selected:
            chrom=site['chrom'];pos=int(site['pos1']);label=f'{chrom}:{pos}'
            old=count_one(a,chrom,pos);new=count_one(b,chrom,pos)
            f=old['fragment_counts'];ref=site['ref'];alt=site['alt']
            comparison={'RNA_fragment_ref':f[ref],'RNA_fragment_alt':f[alt],'RNA_fragment_other':sum(f[c] for c in 'ACGT' if c not in (ref,alt)),'RNA_fragment_discordant':f['conflict'],'RNA_fragment_nonACGT':f['nonACGT']}
            differences={k:{'primary':int(site[k]),'direct':v} for k,v in comparison.items() if int(site[k])!=v}
            assert not differences, f'Independent original-RNA count mismatch at {label}: {differences}'
            oldnames=old['qname_states'];newnames=new['qname_states']
            changed=sorted(n for n in set(oldnames)|set(newnames) if oldnames.get(n)!=newnames.get(n))
            changes=[]
            for n in changed:
                changes.append({'query_name':n,'original_state':oldnames.get(n,'not_qualifying'),'new_state':newnames.get(n,'not_qualifying'),'original_overlapping_records':old['all_overlapping_records_by_name'].get(n,[]),'new_overlapping_records':new['all_overlapping_records_by_name'].get(n,[])})
            row={'chrom':chrom,'pos1':pos,'target_gene':';'.join(within(site)),'ref':ref,'alt':alt,'primary_original_counts_match':True,'changed_qname_states':len(changed),'same_qualifying_qname_states':sum(oldnames[n]==newnames[n] for n in set(oldnames)&set(newnames))}
            for name,obj in [('original',old),('new',new)]:
                for c,v in obj['fragment_counts'].items():row[name+'_'+c]=v
                row[name+'_ACGT']=sum(obj['fragment_counts'][c] for c in 'ACGT')
            row['ACGT_count_vectors_equal']=all(row['original_'+c]==row['new_'+c] for c in 'ACGT')
            row['new_still_ge20_ACGT']=row['new_ACGT']>=20
            results.append(row)
            details[label]={'site':{'chrom':chrom,'pos1':pos,'ref':ref,'alt':alt},'original_qname_states':oldnames,'new_qname_states':newnames,'changed_names':changes,'original_fragment_counts':old['fragment_counts'],'new_fragment_counts':new['fragment_counts'],'original_read_counts':old['read_base_counts'],'new_read_counts':new['read_base_counts'],'original_exclusion_record_counts':old['overlap_record_reasons'],'new_exclusion_record_counts':new['overlap_record_reasons']}
    fields=list(results[0]) if results else ['chrom','pos1','target_gene','ref','alt']
    with (OUT/'comparison.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(results)
    write_json(OUT/'qname-evidence.json',details)
    summary={'status':'PASS','public_loci_in_retained_intervals':len(universe),'public_loci_by_target':dict(collections.Counter(within(x)[0] for x in universe)),'joint_callable_loci_in_retained_intervals':len(selected),'joint_callable_by_target':dict(collections.Counter(within(x)[0] for x in selected)),'loci_with_equal_ACGT_vectors':sum(r['ACGT_count_vectors_equal'] for r in results),'loci_with_changed_qname_states':sum(r['changed_qname_states']>0 for r in results),'total_site_specific_changed_qname_states':sum(r['changed_qname_states'] for r in results),'loci_remaining_at_least20_ACGT_new':sum(r['new_still_ge20_ACGT'] for r in results),'all_independent_original_counts_match_primary':True,'reference_length_checks':reference_checks,'filters':{'MAPQ':30,'BQ':25,'require_proper_pair':True,'exclude_flags':'0xF0C','depth':'uncapped','BAQ':'not recalculated','mate_conflicts':'Qualifying records grouped by query name; >1 A/C/G/T base discards entire name. No UMI-molecule interpretation.'},'limits':['Only preselected jointly callable SNPs within retained BAP1/RASA1 intervals are tested. Other loci cannot be judged from this target BAM. LATS1/LATS2 are outside this six-chromosome public panel.','Retained target BAM contains records intersecting four gene intervals; mates outside intervals were not necessarily retained. A base covering an included SNP intersects its target interval, but full pair geometry elsewhere cannot be reconstructed from this file.','Original and new pipelines differ in STAR version, parameters, sparse-index setting and mapping choices. Differences cannot be assigned to a single factor.','A name not represented at the SNP in the new target BAM may map elsewhere or be unmapped; this is not biological loss of RNA or an allele.','This is a small secondary alignment-dependence check, not a full-panel replication or clinical sample-identity validation.'],'elapsed_seconds':round(time.time()-started,3),'peak_RSS_bytes_macOS':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    summary['primary_coverage_of_target_panel']={'DNA_callable20_loci':sum(int(r['DNA_callable20']) for r in target_coverage),'RNA_callable20_loci':sum(int(r['RNA_callable20']) for r in target_coverage),'maximum_original_RNA_ACGT_depth':max(int(r['RNA_ACGT_depth']) for r in target_coverage),'maximum_DNA_ACGT_depth':max(int(r['DNA_ACGT_depth']) for r in target_coverage),'source':'Existing completed primary counts only; no additional target patient read scan'}
    if not selected:
        summary['status']='NOT_EVALUABLE_NO_QUALIFYING_LOCI'
        summary['all_independent_original_counts_match_primary']=None
        summary['patient_locus_recounts_performed']=0
        summary['interpretation']='No public-panel SNP within the retained target intervals met the frozen joint-callability criterion. No patient per-site recount or QNAME comparison was executed. Zero comparison rows are not evidence of concordance, missing RNA biology, or alignment agreement.'
    write_json(OUT/'summary.json',summary)
    provenance={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':sha(Path(__file__)),'count_complete_sha256':sha(ROOT/'count-complete.json'),'joint_callable_sha256':sha(ROOT/'joint-callable.tsv'),'public_panel_sha256':sha(ROOT/'input/public-common-SNP-panel.tsv'),'original_bam':{'path':str(OLD),'bytes':OLD.stat().st_size,'index':str(OLD_INDEX),'index_sha256':sha(OLD_INDEX),'hash_note':'Previously verified resident original; no redundant full BAM hash pass performed for this bounded check.'},'new_bam':{'path':str(NEW),'bytes':NEW.stat().st_size,'sha256':sha(NEW),'index_sha256':sha(Path(str(NEW)+'.bai'))},'headers':headers,'output_hashes':{p.name:sha(p) for p in [OUT/'comparison.tsv',OUT/'qname-evidence.json',OUT/'summary.json',OUT/'qualification.json']}}
    provenance['all_public_loci_count_source_sha256']=sha(ROOT/'all-public-loci-counts.tsv')
    provenance['output_hashes']['target-panel-coverage.tsv']=sha(OUT/'target-panel-coverage.tsv')
    write_json(OUT/'provenance.json',provenance);print(json.dumps(summary,indent=2))
if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='qualify':print(json.dumps(qualify(),indent=2))
    else:run()
