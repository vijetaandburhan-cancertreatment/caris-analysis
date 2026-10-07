"""Read-only benchmark outcome extraction with exact target coordinates and per-pair evidence."""
from pathlib import Path
import json,csv,gzip,collections,hashlib,datetime,re
import pysam
P=Path(__file__).resolve().parent
design=json.loads((P/'design-manifest.json').read_text());datasets={d['id']:d for d in design['datasets']}
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def table(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def target(r):return r.get('breakpoint1')=='chr22:23290413' and r.get('breakpoint2')=='chr9:130854064' and r.get('#gene1',r.get('gene1'))=='BCR' and r.get('gene2')=='ABL1'
def readnames(rows):return sorted({n for r in rows for n in r.get('read_identifiers','').split(',') if n not in ['', '.']})
allruns=[];allnames=[]
for f in sorted((P/'runs').glob('*/completion.json')):
    meta=json.loads(f.read_text());out=f.parent;d=datasets[meta['dataset']]
    no_input=meta.get('caller_outcome','').startswith('no_chimeric_input')
    accepted=[] if no_input else table(out/'fusions.tsv');discarded=[] if no_input else table(out/'fusions.discarded.tsv');at=[r for r in accepted if target(r)];dt=[r for r in discarded if target(r)]
    junctions=[]
    for line in (out/'STAR.Chimeric.out.junction').read_text().splitlines():
        q=line.split('\t')
        if len(q)<14 or line.startswith('#') or q[0]=='chr_donorA':continue
        junctions.append({'raw_fields':q,'chrom1':q[0],'site1':int(q[1]),'strand1':q[2],'chrom2':q[3],'site2':int(q[4]),'strand2':q[5],'name':q[9],'segment1_cigar':q[11],'segment2_cigar':q[13],'target_junction':{(q[0],int(q[1])),(q[3],int(q[4]))}=={('chr22',23290414),('chr9',130854063)} and q[2]==q[5]})
    tj=sorted({j['name'] for j in junctions if j['target_junction']});anyj=sorted({j['name'] for j in junctions})
    align=collections.defaultdict(list)
    with pysam.AlignmentFile(str(out/'alignments.bam'),'rb') as b:
        for r in b:
            align[r.query_name].append({'flag':r.flag,'mate':1 if r.is_read1 else 2,'reference':r.reference_name,'start0':r.reference_start,'end0':r.reference_end,'cigar':r.cigarstring,'mapq':r.mapping_quality,'mate_reference':r.next_reference_name,'mate_start0':r.next_reference_start,'tlen':r.template_length,'query_length':r.query_length,'is_reverse':r.is_reverse,'is_unmapped':r.is_unmapped,'is_secondary':r.is_secondary,'is_supplementary':r.is_supplementary,'tags':{k:r.get_tag(k) for k in ['NH','HI','NM','nM','AS','SA'] if r.has_tag(k)}})
    assert set(align)==set(d['names']),f'BAM query-name inventory differs from manifest: {out}'
    # Logs may have hundreds of MB of repetitive absent-contig warnings; stream, never expand on disk.
    cache=P/'analysis-cache'/f'{meta["id"]}.warnings.json';cache.parent.mkdir(exist_ok=True)
    log_hash=sha(out/'Arriba.log.gz')
    if cache.exists():
        cached=json.loads(cache.read_text());assert cached['log_sha256']==log_hash
        warnings=cached['warning_categories'];malformed=cached['malformed_warning_lines']
    else:
        warnings=collections.Counter();malformed=[]
        needles={'warning_tokens':b'WARNING:', 'unknown_gene_or_malformed_range_tokens':b'WARNING: unknown gene or malformed range:', 'default_mate_gap_warning_tokens':b'WARNING: not enough chimeric reads to estimate mate gap distribution'}
        tail=b''
        with gzip.open(out/'Arriba.log.gz','rb') as log:
            while True:
                chunk=log.read(65536);data=tail+chunk;cut=max(0,len(data)-512) if chunk else len(data)
                remain=data[cut:]
                for key,needle in needles.items():warnings[key]+=data.count(needle)-remain.count(needle)
                for hit in re.finditer(rb'WARNING: [0-9]+ SAM records were malformed and ignored',data):
                    if hit.start()<cut:malformed.append(hit.group().decode())
                tail=remain
                if not chunk:break
        cache.write_text(json.dumps({'method':'bounded64KiB chunks; literal token occurrences, not unique events','log_sha256':log_hash,'warning_categories':dict(warnings),'malformed_warning_lines':malformed},indent=2))
    for name in d['names']:
        in_at=name in readnames(at);in_dt=name in readnames(dt)
        if in_at:state='exact_target_accepted_support'
        elif name in tj:state='STAR_exact_target_junction_not_in_accepted_target_support'
        elif name in anyj:state='other_STAR_chimeric_junction_not_exact_target'
        elif all(x['is_unmapped'] for x in align[name]):state='no_mapped_record'
        else:state='mapped_without_exact_target_STAR_chimeric_junction'
        allnames.append({'run':meta['id'],'dataset':d['id'],'name':name,'source_fusion_pair':name in d['source_fusion_names'],'STAR_exact_target_junction':name in tj,'Arriba_accepted_exact_target_support':in_at,'Arriba_discarded_exact_target_support':in_dt,'state':state,'target_discarded_filters':[r['filters'] for r in dt if name in r.get('read_identifiers','').split(',')],'alignments':align[name],'chimeric_junctions':[j for j in junctions if j['name']==name]})
    row={'caller_outcome':meta['caller_outcome'],'run':meta['id'],'dataset':d['id'],'whole_reference':meta['whole_reference'],'sparseD':meta['sparseD'],'known_fusion_recovery':meta['known_fusion_recovery'],'input_pairs':d['input_pairs'],'selected_source_pairs':len(d['source_fusion_names']),'target_accepted':bool(at),'target_confidence':[r['confidence'] for r in at],'accepted_support_names':readnames(at),'accepted_support_count_names':len(readnames(at)),'STAR_target_junction_names':tj,'STAR_target_junction_count_names':len(tj),'target_discarded_filters':[r['filters'] for r in dt],'accepted_target_rows':at,'discarded_target_rows':dt,'other_accepted_rows':[r for r in accepted if not target(r)],'all_chimeric_names':anyj,'malformed_warning_lines':malformed,'warning_categories':dict(warnings),'STAR_resources':meta['STAR_resources'],'Arriba_resources':meta['Arriba_resources'],'completion_sha256':sha(f)}
    allruns.append(row)
by={(r['dataset'],r['sparseD']):r for r in allruns if not r['whole_reference'] and r['known_fusion_recovery']}
comparisons=[]
for name in datasets:
    if (name,1) not in by or (name,8) not in by:continue
    a,b=by[(name,1)],by[(name,8)]
    comparisons.append({'dataset':name,'accepted_D1':a['target_accepted'],'accepted_D8':b['target_accepted'],'STAR_names_only_D1':sorted(set(a['STAR_target_junction_names'])-set(b['STAR_target_junction_names'])),'STAR_names_only_D8':sorted(set(b['STAR_target_junction_names'])-set(a['STAR_target_junction_names'])),'Arriba_names_only_D1':sorted(set(a['accepted_support_names'])-set(b['accepted_support_names'])),'Arriba_names_only_D8':sorted(set(b['accepted_support_names'])-set(a['accepted_support_names']))})
out={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'design_manifest_sha256':sha(P/'design-manifest.json'),'runs_completed':len(allruns),'expected_runs':design['resources']['expected_max_runs'],'runs':allruns,'paired_compact_comparisons':comparisons,'interpretation_limits':[x.replace('only three prespecified cases','only four prespecified cases') for x in design['limitations']]+['The fixed32-parent synthetic background cannot calibrate caller confidence, e-values, recovery thresholds, RNA abundance or detection limits for the23million-pair patient library.','Per-name states describe observed evidence, not a causal explanation of every rejection. Arriba filter labels apply to candidate hypotheses; no per-read biological false-negative rate is inferred.','STAR junction coordinates and Arriba gene breakpoints use different boundary conventions; exact conversions are declared in the parser.']}
(P/'results.json').write_text(json.dumps(out,indent=2));(P/'per-pair-evidence.json').write_text(json.dumps(allnames,indent=2))
fields=['run','caller_outcome','dataset','whole_reference','sparseD','known_fusion_recovery','input_pairs','selected_source_pairs','target_accepted','target_confidence','accepted_support_count_names','accepted_support_names','STAR_target_junction_count_names','STAR_target_junction_names','target_discarded_filters']
with (P/'results.tsv').open('w') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader()
    for r in allruns:w.writerow({k:','.join(map(str,r[k])) if isinstance(r[k],list) else r[k] for k in fields})
print(json.dumps({'runs_completed':len(allruns),'expected':design['resources']['expected_max_runs'],'target_accepted_cases':sum(r['target_accepted'] for r in allruns),'compared_datasets':len(comparisons)},indent=2))
