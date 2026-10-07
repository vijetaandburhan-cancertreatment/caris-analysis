"""Separate read-only outcome parser for the prespecified modest public fusion study."""
from pathlib import Path
import json, hashlib, gzip, re, collections, datetime, argparse
import pysam
P=Path(__file__).resolve().parent
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rc(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def fq(p):
    with gzip.open(p,'rt') as f:
        while h:=f.readline():
            s=f.readline().strip();plus=f.readline();q=f.readline().strip()
            assert h.startswith('@') and plus.startswith('+') and len(s)==len(q)
            yield re.sub(r'/[12]$','',h[1:].split()[0]),s,q
def table(p):
    lines=p.read_text().splitlines();keys=lines[0].split('\t');out=[]
    for line in lines[1:]:
        fields=line.split('\t');assert len(fields)==len(keys)
        out.append(dict(zip(keys,fields)))
    return out
def exact_target(row):
    # Reversed biological event would not equal this directed fusion.
    return row['#gene1']=='BCR' and row['gene2']=='ABL1' and row['breakpoint1']=='chr22:23290413' and row['breakpoint2']=='chr9:130854064'
def names(rows):return set(n for x in rows for n in x['read_identifiers'].split(',') if n not in ['','.'])
def diagnostics(path):
    # Compact-reference logs can contain enormous progress/warning lines.
    # Scan bounded chunks; preserve short diagnostics, never whole unbounded lines.
    patterns={'malformed':re.compile(rb'[0-9]+ SAM records were malformed'),
              'error':re.compile(rb'ERROR: [^\r\n]{1,250}')}
    found={k:[] for k in patterns};seen={k:set() for k in patterns};carry=b'';offset=0
    with gzip.open(path,'rb') as f:
        while block:=f.read(65536):
            window=carry+block;start=offset-len(carry)
            for kind,pattern in patterns.items():
                for m in pattern.finditer(window):
                    position=start+m.start()
                    if position not in seen[kind]:found[kind].append(m.group().decode(errors='replace'));seen[kind].add(position)
            offset+=len(block);carry=window[-512:]
    return found
def bam_audit(path,expected_lengths):
    records=collections.defaultdict(list);touch=collections.defaultdict(lambda:collections.defaultdict(lambda:collections.defaultdict(list)))
    with pysam.AlignmentFile(str(path),'rb') as bam:
        for r in bam:
            mate=1 if r.is_read1 else 2 if r.is_read2 else 0
            assert mate in [1,2],(r.query_name,r.flag)
            summary={'mate':mate,'flag':r.flag,'reference':r.reference_name,'start0':r.reference_start,'end0':r.reference_end,'cigar':r.cigarstring,'MAPQ':r.mapping_quality,'NH':r.get_tag('NH') if r.has_tag('NH') else None,'HI':r.get_tag('HI') if r.has_tag('HI') else None}
            records[r.query_name].append(summary)
            if r.is_unmapped or r.reference_name not in ['chr22','chr9']:continue
            role='BCR' if r.reference_name=='chr22' else 'ABL1';pos0=23290412 if role=='BCR' else 130854063
            lead_hard=r.cigartuples[0][1] if r.cigartuples and r.cigartuples[0][0]==5 else 0
            for q,g in r.get_aligned_pairs(matches_only=True):
                if g==pos0:
                    coord=lead_hard+q
                    if r.is_reverse:coord=expected_lengths[r.query_name][mate-1]-1-coord
                    assert 0<=coord<expected_lengths[r.query_name][mate-1]
                    touch[r.query_name][mate][role].append({'original_read_base0':coord,'flag':r.flag,'MAPQ':r.mapping_quality,'HI':summary['HI']})
    adjacent={}
    for name,ms in touch.items():
        events=[]
        for mate,roles in ms.items():
            for a in roles.get('BCR',[]):
                for b in roles.get('ABL1',[]):
                    same_orientation=bool(a['flag']&16)==bool(b['flag']&16)
                    consistent_order=b['original_read_base0']-a['original_read_base0']==(-1 if a['flag']&16 else 1)
                    same_hit=a['HI']==b['HI'] or a['HI'] is None or b['HI'] is None
                    if same_orientation and consistent_order and same_hit:
                        events.append({'mate':mate,'BCR':a,'ABL1':b})
        if events:adjacent[name]=events
    return dict(records),adjacent

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--partial',action='store_true');ap.add_argument('--reuse-reviewed-partial',action='store_true');args=ap.parse_args()
    design=json.loads((P/'design-manifest.json').read_text());datasets={x['id']:x for x in design['datasets']};manifest=sha(P/'design-manifest.json')
    assert json.loads((P/'independent-design-review.json').read_text())['manifest_sha256']==manifest
    fusion=''.join((P/'truth-fusion-transcript.fa').read_text().splitlines()[1:]);j=design['junction']['fusion_transcript_boundary0']
    input_checks={}
    for name,d in datasets.items():
        aa,bb=[list(fq(f['path'])) for f in d['files']];assert len(aa)==len(bb)==d['input_pairs']
        pairs={a[0]:(a,b) for a,b in zip(aa,bb)};assert all(a[0]==b[0] for a,b in pairs.values())
        exact={}
        for anchor in [12,15,20,25,30]:
            marker=fusion[j-anchor:j+anchor];patterns=[marker,rc(marker)]
            for qmin in [0,20,30]:
                support=[]
                for n,pp in pairs.items():
                    found=False
                    for _,s,q in pp:
                        for pattern in patterns:
                            at=s.find(pattern)
                            while at>=0:
                                if min(ord(x)-33 for x in q[at:at+len(pattern)])>=qmin:found=True
                                at=s.find(pattern,at+1)
                    if found:support.append(n)
                exact[f'{anchor}+{anchor}_Q{qmin}']=sorted(support)
        input_checks[name]={'pair_lengths':{n:[len(p[0][1]),len(p[1][1])] for n,p in pairs.items()},'unique_pair_names':len(pairs),'unique_full_sequence_pairs':len(set((p[0][1],p[1][1]) for p in pairs.values())),'exact_reference_junction_marker_support':exact}
    runs=[];cache={};binding=None;reused=[]
    if args.reuse_reviewed_partial:
        assert not args.partial,'Do not overwrite the bound partial while using it.'
        binding=json.loads((P/'independent-partial-cache-binding.json').read_text())
        assert sha(P/'independent-outcome-partial.json')==binding['partial_sha256']
        previous=json.loads((P/'independent-outcome-partial.json').read_text())
        assert previous['manifest_sha256']==binding['manifest_sha256']==manifest
        assert previous['audit_script_sha256']==binding['producer_audit_script_sha256']
        cache={r['id']:r for r in previous['runs']}
    for path in sorted((P/'runs').glob('*/completion.json')):
        meta=json.loads(path.read_text());out=path.parent;d=datasets[meta['dataset']]
        assert meta['status']=='complete' and meta['manifest_sha256']==manifest
        assert meta['STAR_resources']['exit_codes']==[0]
        no_input=meta.get('caller_outcome','').startswith('no_chimeric_input;')
        assert meta['Arriba_resources']['exit_codes']==([1,0] if no_input else [0,0])
        for filename,item in meta['outputs'].items():assert sha(out/filename)==item['sha256']
        if meta['id'] in cache:
            assert sha(path)==binding['completion_file_sha256'][meta['id']]
            runs.append(cache[meta['id']]);reused.append(meta['id']);continue
        if no_input:
            assert not (out/'fusions.tsv').exists() and not (out/'fusions.discarded.tsv').exists()
            accepted=[];discarded=[]
        else:accepted=table(out/'fusions.tsv');discarded=table(out/'fusions.discarded.tsv')
        at=[x for x in accepted if exact_target(x)];dt=[x for x in discarded if exact_target(x)]
        an,dn=names(at),names(dt);assert (an|dn)<=set(d['names'])
        chim=[]
        for line in (out/'STAR.Chimeric.out.junction').read_text().splitlines():
            f=line.split('\t')
            if len(f)<14 or not f[1].lstrip('-').isdigit():continue
            # Convert STAR intronic boundaries explicitly in each directed orientation.
            forward=(f[0],int(f[1])-1,f[2],f[3],int(f[4])+1,f[5])==('chr22',23290413,'+','chr9',130854064,'+')
            reverse=(f[0],int(f[1])+1,f[2],f[3],int(f[4])-1,f[5])==('chr9',130854064,'-','chr22',23290413,'-')
            chim.append({'name':f[9],'target':forward or reverse,'chromosomes':[f[0],f[3]],'boundaries_STAR1':[int(f[1]),int(f[4])],'strands':[f[2],f[5]],'segment_starts_STAR1':[int(f[10]),int(f[12])],'CIGARs':[f[11],f[13]],'repeat_lengths':[int(f[7]),int(f[8])]})
        tn={x['name'] for x in chim if x['target']};assert tn<=set(d['names'])
        records,adjacent=bam_audit(out/'alignments.bam',input_checks[d['id']]['pair_lengths'])
        assert set(records)==set(d['names'])
        no_input_check=None
        if no_input:
            assert not chim and int(meta['STAR_final']['Number of chimeric reads'])==0
            counts=collections.Counter();allnames=set()
            pysam.samtools.quickcheck(str(out/'alignments.bam'))
            with pysam.AlignmentFile(str(out/'alignments.bam'),'rb') as bam:
                for r in bam:
                    allnames.add(r.query_name);counts['records']+=1;counts['supplementary']+=r.is_supplementary;counts['secondary']+=r.is_secondary;counts['unmapped']+=r.is_unmapped;counts['improper']+=not r.is_proper_pair;counts['SA']+=r.has_tag('SA');counts['different_mate_chromosome']+=not r.is_unmapped and not r.mate_is_unmapped and r.reference_id!=r.next_reference_id
            assert allnames==set(d['names']) and counts['supplementary']==counts['SA']==counts['different_mate_chromosome']==0
            no_input_check={'status':'explicit no-input terminal state verified, no call tables generated','distinct_query_names':len(allnames),'BAM_counts':dict(counts),'STAR_chimeric_rows':len(chim),'BAM_quickcheck':'PASS'}
        for r in at:
            total=sum(int(r[k]) for k in ['split_reads1','split_reads2','discordant_mates'])
            assert total==len(names([r])),('caller counts vs unique-name disagreement',meta['id'])
        cmd=meta['Arriba_command'];assert ('-k' in cmd)==meta['known_fusion_recovery'] and '-t' in cmd and '-f' not in cmd
        assert meta['STAR_command'][meta['STAR_command'].index('--runThreadN')+1]=='2'
        if not meta['known_fusion_recovery']:assert meta['dataset']==design['no_prior_comparison']['dataset'] and meta['sparseD']==1 and not meta['whole_reference']
        # Parse warnings separately; do not collapse the mixed-unit Arriba counter to input reads.
        log_diagnostics=diagnostics(out/'Arriba.log.gz');malformed=log_diagnostics['malformed'];other_warnings={}
        no_input_error=[x for x in log_diagnostics['error'] if x.startswith('ERROR: no split reads or discordant mates found')]
        other_errors=[x for x in log_diagnostics['error'] if not x.startswith('ERROR: no split reads or discordant mates found')]
        assert not other_errors and bool(no_input_error)==no_input
        runs.append({'id':meta['id'],'dataset':d['id'],'whole_reference':meta['whole_reference'],'sparseD':meta['sparseD'],'known_fusion_recovery':meta['known_fusion_recovery'],'caller_outcome':meta.get('caller_outcome'),'no_chimeric_input_independent_validation':no_input_check,'no_input_error_lines':no_input_error,
          'input_pairs':len(records),'accepted_row_count':len(accepted),'discarded_row_count':len(discarded),'exact_target_accepted_rows':at,'exact_target_discarded_rows':dt,'accepted_target_names':sorted(an),'discarded_target_names':sorted(dn),'STAR_exact_target_names':sorted(tn),'separate_BAM_adjacent_junction_names':sorted(adjacent),'BAM_adjacent_junction_evidence':adjacent,'STAR_target_not_caller_accepted':sorted(tn-an),'caller_accepted_not_STAR_exact':sorted(an-tn),'chimeric_rows':chim,'BAM_records_by_name':records,'malformed_warning_lines':malformed,'other_warning_lines':dict(other_warnings)})
    if not args.partial:assert len(runs)==design['resources']['expected_max_runs']==45
    by={(r['dataset'],r['sparseD']):r for r in runs if not r['whole_reference'] and r['known_fusion_recovery']};comparisons=[]
    for dataset in datasets:
        if (dataset,1) not in by or (dataset,8) not in by:continue
        a,b=by[dataset,1],by[dataset,8]
        comparisons.append({'dataset':dataset,'D1_accepted_names':a['accepted_target_names'],'D8_accepted_names':b['accepted_target_names'],'D1_only':sorted(set(a['accepted_target_names'])-set(b['accepted_target_names'])),'D8_only':sorted(set(b['accepted_target_names'])-set(a['accepted_target_names'])),'D1_STAR_only':sorted(set(a['STAR_exact_target_names'])-set(b['STAR_exact_target_names'])),'D8_STAR_only':sorted(set(b['STAR_exact_target_names'])-set(a['STAR_exact_target_names']))})
    summary_comparison={'status':'owner results not yet available'}
    if (P/'results.json').exists():
        owner=json.loads((P/'results.json').read_text());ob={x['run']:x for x in owner['runs']};compared=0
        for r in runs:
            if r['id'] not in ob:continue
            x=ob[r['id']]
            assert r['accepted_target_names']==x['accepted_support_names']
            assert r['STAR_exact_target_names']==x['STAR_target_junction_names']
            assert bool(r['exact_target_accepted_rows'])==x['target_accepted']
            assert [x['filters'] for x in r['exact_target_discarded_rows']]==x['target_discarded_filters']
            compared+=1
        summary_comparison={'status':'PASS','runs_compared':compared}
    result={'status':'PARTIAL' if args.partial else 'PASS','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'manifest_sha256':manifest,'audit_script_sha256':sha(__file__),'completed_runs_independently_audited':len(runs),'planned_runs':45,'independent_input_marker_counts':input_checks,'runs':runs,'paired_compact_comparisons':comparisons,'comparison_to_primary_extractor':summary_comparison,
     'methods':'Separate TSV parsing and explicit STAR-boundary conversion; direct BAM aligned-pair reconstruction of adjacent retained BCR/ABL1 bases in original mate coordinates; exact raw-input marker+base-quality recount; unique-query-name collapse; independent output hash, exit code, input inventory and command checks.',
     'limitations':['Public synthetic fixture only, not patient validation or general sensitivity/LoD.','Exact-reference marker absence is not failure of a known true fusion: mismatches, trimming and base-quality cutoffs can eliminate support.','Direct BAM adjacency is an additional local alignment representation check; it need not cover every caller-supported geometry or paired-end merged representation.','Matched subsets overlap; distinct query names or sequence families do not establish independent biological molecules.','Counter warnings have mixed query-name/HI-group and alignment-record units; no percentage of lost input pairs can be inferred.','Compact D1/D8 comparisons lack genome-wide competitors, and full-D8 controls do not establish full-D1 equality.']}
    result['verified_partial_reuse']={'reused_runs':reused,'binding_sha256':sha(P/'independent-partial-cache-binding.json') if binding else None,'producer_script_sha256':binding['producer_audit_script_sha256'] if binding else None,'validation':'Every original output file hash and exact completion-manifest hash rechecked before reusing already completed independent tests.'}
    destination=P/('independent-outcome-partial.json' if args.partial else 'independent-outcome-audit.json');destination.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'runs':len(runs),'owner_comparison':summary_comparison}))
if __name__=='__main__':main()
