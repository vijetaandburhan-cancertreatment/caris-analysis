"""Fast complete scientific-claim check against hash-verified resident immutable calls/BAMs."""
from pathlib import Path
import importlib.util,json,hashlib,datetime,pysam
P=Path(__file__).resolve().parent
R=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct5-fusion-robustness-v1')
s=importlib.util.spec_from_file_location('independent',P/'independent_outcome_audit.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
d=json.loads((P/'design-manifest.json').read_text());ds={x['id']:x for x in d['datasets']};primary=json.loads((P/'results.json').read_text());assert primary['runs_completed']==45
lengths={}
for name,x in ds.items():
    aa,bb=[list(m.fq(R/'inputs'/Path(f['path']).name)) for f in x['files']]
    lengths[name]={a[0]:[len(a[1]),len(b[1])] for a,b in zip(aa,bb)}
    assert len(lengths[name])==x['input_pairs']
out=[]
for pr in primary['runs']:
    p=R/'runs'/pr['run'];meta=json.loads((p/'completion.json').read_text());x=ds[meta['dataset']]
    noinput=meta['caller_outcome'].startswith('no_chimeric_input')
    assert meta['STAR_resources']['exit_codes']==[0] and meta['Arriba_resources']['exit_codes']==([1,0] if noinput else [0,0])
    nameset=set(x['names']);records,adj=m.bam_audit(p/'alignments.bam',lengths[x['id']]);assert set(records)==nameset
    ac=[] if noinput else m.table(p/'fusions.tsv');dc=[] if noinput else m.table(p/'fusions.discarded.tsv')
    at=[r for r in ac if m.exact_target(r)];dt=[r for r in dc if m.exact_target(r)];an=sorted(m.names(at));tn=set()
    for line in (p/'STAR.Chimeric.out.junction').read_text().splitlines():
        q=line.split('\t')
        if len(q)<14 or not q[1].isdigit():continue
        direct=(q[0],int(q[1])-1,q[2],q[3],int(q[4])+1,q[5])==('chr22',23290413,'+','chr9',130854064,'+')
        reverse=(q[0],int(q[1])+1,q[2],q[3],int(q[4])-1,q[5])==('chr9',130854064,'-','chr22',23290413,'-')
        if direct or reverse:tn.add(q[9])
    assert an==pr['accepted_support_names'] and sorted(tn)==pr['STAR_target_junction_names']
    assert bool(at)==pr['target_accepted'] and [r['filters'] for r in dt]==pr['target_discarded_filters']
    assert len(an)==pr['accepted_support_count_names'] and len(tn)==pr['STAR_target_junction_count_names']
    assert sorted(adj)==sorted(tn),'Separate original-mate BAM reconstruction differs from exact STAR names'
    for row in at:assert sum(int(row[k]) for k in ['split_reads1','split_reads2','discordant_mates'])==len(m.names([row]))
    verified_files=[]
    for name in ['alignments.bam','STAR.Chimeric.out.junction','STAR.Log.final.out','fusions.tsv','fusions.discarded.tsv']:
        if name in meta['outputs']:
            assert m.sha(p/name)==meta['outputs'][name]['sha256'];verified_files.append(name)
    if noinput:
        assert not (p/'fusions.tsv').exists() and not (p/'fusions.discarded.tsv').exists() and not tn
        z=m.diagnostics(p/'Arriba.log.gz');assert len(z['error'])==1 and z['error'][0].startswith('ERROR: no split reads or discordant mates found')
        with pysam.AlignmentFile(str(p/'alignments.bam'),'rb') as bam:
            for b in bam:assert not b.is_supplementary and not b.has_tag('SA') and (b.is_unmapped or b.mate_is_unmapped or b.reference_id==b.next_reference_id)
    out.append({'run':pr['run'],'caller_outcome':meta['caller_outcome'],'input_pairs':len(records),'accepted_names':an,'STAR_exact_names':sorted(tn),'independent_BAM_adjacent_junction_names':sorted(adj),'discarded_candidate_filters':[r['filters'] for r in dt],'hash_verified_material_files':verified_files})
assert len(out)==45
by={r['run']:r for r in out};changed=[]
for name in ds:
    a,b=by['compactD1__'+name],by['compactD8__'+name]
    assert bool(a['accepted_names'])==bool(b['accepted_names'])
    if a['accepted_names']!=b['accepted_names']:
        assert set(a['accepted_names'])-set(b['accepted_names'])=={'BCR-ABL1-28'} and not set(b['accepted_names'])-set(a['accepted_names']);changed.append(name)
assert len(changed)==4
assert [len(by[k+'__support_all']['accepted_names']) for k in ['compactD1','compactD8','fullD8']]==[10,9,10]
assert set(by['compactD1__support_all']['accepted_names'])-set(by['fullD8__support_all']['accepted_names'])=={'BCR-ABL1-28'}
assert set(by['fullD8__support_all']['accepted_names'])-set(by['compactD1__support_all']['accepted_names'])=={'BCR-ABL1-18'}
assert sum(r['caller_outcome'].startswith('no_chimeric_input') for r in out)==7
result={'status':'PASS','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'manifest_sha256':m.sha(P/'design-manifest.json'),'primary_results_sha256':m.sha(P/'results.json'),'primary_summary_sha256':m.sha(P/'summary.json'),'reviewed_findings_sha256':m.sha(P/'findings.txt'),'script_sha256':m.sha(__file__),'source_root':str(R),'runs_checked':len(out),'runs':out,'support_identity_difference_conditions':changed,'conclusion':'All45 exact target decisions, accepted query-name counts/identities, candidate-level discarded labels, STAR junction identities and separate BAM junction reconstruction agree.38 ordinary caller-table runs and7 narrowly verified no-input outcomes. Scientific claims in current corrected findings and summary pass independent review.','qualifications':['20 condition labels comprise19 distinct inputs; correlated synthetic data, not clinical sensitivity/specificity or limit-of-detection.','Paired trimming can restore long merged arms; no isolated12nt effective-anchor experiment.','Full-reference D8 versus compact changes reference competition; no full D1/D8 equivalence test.','Comprehensive audit separately rechecks compressed logs/output hashes; this claims audit checks all material call/BAM/junction hashes and exact no-input logs, not every long benign warning.','No new patient analysis or treatment inference.']}
(P/'independent-final-claims-review.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':'PASS','runs_checked':45,'claim_review':'approved','size_bytes':(P/'independent-final-claims-review.json').stat().st_size}))
