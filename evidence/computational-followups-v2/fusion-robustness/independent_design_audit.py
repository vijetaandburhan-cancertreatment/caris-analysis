"""Independent public-input audit: seeded exhaustive <=3-Hamming placements, no owner code imports."""
from pathlib import Path
import gzip, json, hashlib, re, datetime
import pysam

P=Path(__file__).resolve().parent
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rc(s): return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def read_fastq(p):
    rows=[]
    with gzip.open(p,'rt') as f:
        while True:
            h=f.readline()
            if not h: break
            s=f.readline().rstrip('\n\r'); plus=f.readline();q=f.readline().rstrip('\n\r')
            assert h.startswith('@') and plus.startswith('+') and s and len(s)==len(q)
            name=re.sub(r'/[12]$','',h[1:].split()[0])
            rows.append((name,s,q))
    return rows
def pairs(a,b):
    aa,bb=read_fastq(a),read_fastq(b); assert len(aa)==len(bb)
    out={}
    for x,y in zip(aa,bb):
        assert x[0]==y[0] and x[0] not in out
        out[x[0]]=(x[1:],y[1:])
    assert len(set((x[0][0],x[1][0]) for x in out.values()))==len(out)
    return out
def placements(sequence,reference):
    # Pigeonhole guarantee: any ungapped alignment with <=3 mismatches has
    # an exact one of four disjoint seed blocks. All candidate offsets are scored.
    out=[];n=len(sequence)
    for ori,query in [('+',sequence),('-',rc(sequence))]:
        candidates=set()
        for i in range(4):
            a,b=i*n//4,(i+1)*n//4; seed=query[a:b]
            at=reference.find(seed)
            while at>=0:
                start=at-a
                if 0<=start<=len(reference)-n: candidates.add(start)
                at=reference.find(seed,at+1)
        for start in candidates:
            mm=sum(a!=b for a,b in zip(query,reference[start:start+n]))
            if mm<=3:out.append((mm,ori,start))
    if not out:return []
    best=min(x[0] for x in out)
    return sorted(x for x in out if x[0]==best)

d=json.loads((P/'design-manifest.json').read_text());freeze=json.loads((P/'design-freeze.json').read_text())
manifest_hash=sha(P/'design-manifest.json');assert manifest_hash==freeze['manifest_sha256']
assert manifest_hash=='4abf92bfee2d29e0b68f60513f9a5083a99bacb1fb9ec84c4b10665a9808ca09'
assert not (P/'execution-start.json').exists(),'This audit must precede execution.'
fa=pysam.FastaFile(d['source_provenance']['reference_fasta']['path'])
gtf=Path(d['source_provenance']['compact_gtf']['path'])
exons={g:[] for g in d['parent_transcripts']}
for line in gtf.read_text().splitlines():
    if line.startswith('#'):continue
    f=line.split('\t')
    if len(f)!=9 or f[2]!='exon':continue
    attrs=dict(re.findall(r'(\S+) "([^"]*)";',f[8]))
    for g,p in d['parent_transcripts'].items():
        if attrs.get('transcript_id')==p['tx']:
            assert f[6]=='+' and f[0]==p['chrom']
            exons[g].append((int(f[3]),int(f[4])))
parents={};cuts={}
for g,p in d['parent_transcripts'].items():
    ex=sorted(set(exons[g])); coords=[];seqs=[]
    # Independent explicit transcript genomic coordinate map.
    for a,b in ex:
        coords.extend(range(a,b+1));seqs.append(fa.fetch(p['chrom'],a-1,b).upper())
    parents[g]=''.join(seqs)
    cuts[g]=coords.index(p['breakpoint1'])+(1 if g=='BCR' else 0)
    assert cuts[g]==p['cut0']
    assert hashlib.sha256(parents[g].encode()).hexdigest()==p['sequence_sha256']
    assert [[a-1,b] for a,b in ex]==p['exons0']
fusion=parents['BCR'][:cuts['BCR']]+parents['ABL1'][cuts['ABL1']:];j=cuts['BCR']
assert hashlib.sha256(fusion.encode()).hexdigest()==d['junction']['fusion_transcript_sha256']
assert j==d['junction']['fusion_transcript_boundary0']
source=pairs(d['source_provenance']['source_R1']['path'],d['source_provenance']['source_R2']['path'])
assert len(source)==40
truth={};pool=[];placements_checked=0
owner_truth={x['name']:x for x in d['truth_reads']}
for name,rs in source.items():
    tr=[]
    for mi,(seq,quality) in enumerate(rs):
        hits=placements(seq,fusion);own=owner_truth[name]['mates'][mi]
        good=False
        if hits:
            assert len(hits)==own['best_locations'] and hits[0][0]==own['best_mismatches']
            assert hits[0][1:]==(own['orientation'],own['start_on_fusion0'])
            placements_checked+=1
        else: assert own['best_mismatches']>3
        if len(hits)==1:
            mm,ori,start=hits[0];n=len(seq)
            crossing=start<j<start+n
            left,right=(j-start,start+n-j) if crossing else (None,None)
            good=crossing and min(left,right)>=10
            tr.append({'orientation':ori,'start':start,'mismatches':mm,'crossing':crossing,
                       'left':left,'right':right,'cut':(left if ori=='+' else right) if crossing else None,'eligible':good})
        else:tr.append({'eligible':False,'crossing':False})
        assert good==own['eligible_crossing_mate']
    truth[name]=tr
    if any(x['eligible'] for x in tr):pool.append(name)
assert pool==d['eligible_source_names'] and len(pool)==13
def minor(name):return min(min(x['left'],x['right']) for x in truth[name] if x['eligible'])
strong=sorted((n for n in pool if minor(n)>=40),key=lambda n:(-minor(n),n))[:4]
assert strong==d['strong_anchor_names']
assert sorted(pool,key=lambda n:(minor(n),n))[0]==d['shortest_anchor_name']
background={}
for row in d['parent_fragments']:
    a,b=row['fragment_start0'],row['fragment_end0'];par=parents[row['gene']]
    assert b-a==row['fragment_length'] and 0<=a<b<=len(par)
    qa,qb=source[row['quality_source_name']]
    background[row['name']]=((par[a:a+150],qa[1]),(rc(par[b-150:b]),qb[1]))
assert len(background)==32
markers={str(n):fusion[j-n:j+n] for n in [15,20,25]}
negative_marker_names={n:[name for name,p in background.items() if any(mark in seq or rc(mark) in seq for seq,q in p)] for n,mark in markers.items()}
assert all(not x for x in negative_marker_names.values())
checks=[]
for ds in d['datasets']:
    actual=pairs(*(x['path'] for x in ds['files']))
    assert len(actual)==ds['input_pairs'] and list(actual)==ds['names']
    for f in ds['files']:
        assert sha(f['path'])==f['sha256']
        assert hashlib.sha256(gzip.decompress(Path(f['path']).read_bytes())).hexdigest()==f['uncompressed_sha256']
    selected=ds['source_fusion_names'];expected=dict(source) if ds['id']=='fixture_all40' else dict(background)
    if ds['id']!='fixture_all40':
        for name in selected:
            assert name in pool
            expected[name]=source[name]
    if 'trim_records' in ds:
        cap=ds['anchor_cap'];mrecords={(r['name'],r['mate']):r for r in ds['trim_records']}
        assert selected==strong
        for name in selected:
            changed=[]
            for mi,((seq,q),tr) in enumerate(zip(source[name],truth[name]),1):
                a,b=0,len(seq)
                if cap is not None and tr['crossing']:
                    cut=tr['cut']
                    if cut<=len(seq)-cut:a=cut-cap
                    else:b=cut+cap
                    assert min(cut-a,b-cut)==cap
                m=mrecords[name,mi]
                assert [a,b,b-a]==[m['original_start0_retained'],m['original_end0_retained'],m['retained_length']]
                changed.append((seq[a:b],q[a:b]))
            expected[name]=tuple(changed)
    assert actual==expected,ds['id']
    if ds['id'].startswith('support_') and ds['id']!='support_all':
        seed=ds['seed'];n=ds['source_support_level']
        ranked=sorted(pool,key=lambda name:hashlib.sha256(f'fusion-robustness-v1:seed{seed}:{name}'.encode()).hexdigest())
        assert selected==sorted(ranked[:n]) and len(selected)==n
    checks.append({'dataset':ds['id'],'pairs':len(actual),'unique_names':len(actual),'unique_full_sequence_pairs':len(set((r[0][0],r[1][0]) for r in actual.values())),'literal_sequence_and_quality_comparison':'PASS'})
assert len(checks)==20
assert d['compact_matrix']==[{'dataset':x['id'],'sparseD':s} for x in d['datasets'] for s in [1,8]]
assert d['full_reference_D8_transfer_checks']==['parental_zero','support_all','strong_original','strong_trim_anchor12']
assert d['resources']['expected_max_runs']==45
assert d['no_prior_comparison']['omit_arriba_flags']==['-k'] and d['no_prior_comparison']['keep_tags_flag_t']
review={'status':'approved','manifest_sha256':manifest_hash,'reviewed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
 'reviewer':'independent subagent caris_variants_inventory','audit_script_sha256':sha(__file__),
 'scope':'Pre-execution public-data study only; no new patient analysis or clinical validity claim.',
 'independent_methods':['Separate stdlib FASTQ parser and paired EOF/name/sequence/quality checks.','Independent GTF parsing and explicit genomic-base-to-transcript coordinate map.','Four-seed pigeonhole enumeration of every <=3-mismatch ungapped placement, separate from owner numpy window distances.','Literal independent reconstruction of all parental fragments and trim intervals.'],
 'reference_cuts':cuts,'fusion_sha256':hashlib.sha256(fusion.encode()).hexdigest(),'source_pairs':len(source),'eligible_pairs':len(pool),'mates_with_at_most_3_mismatches_checked':placements_checked,
 'eligible_names':pool,'strong_names':strong,'shortest_anchor_name':d['shortest_anchor_name'],'parental_negative_marker_names':negative_marker_names,'datasets':checks,
 'runner_review':'Exact manifest/input hash gate; same STAR COMMON and stock Arriba filters; -k removed only in the prespecified no-prior run while -t retained; sequential runs, resource guards and preserved logs.',
 'limitations':['Subsets overlap and are not independent biological replicates.','Thirteen eligible pairs define a selected public fixture pool, not every possible true fusion read.','No duplicate pair-name or full-sequence inflation within any dataset; source names recur across planned comparisons.','Compact reference omits competitors; only four small full-D8 comparisons, with no full-D1 counterpart.','Trimming changes length and mate geometry as well as anchor size; those effects cannot be isolated.','Thirty-two synthetic normal-parent fragments are an execution negative control, not a transcriptome-wide specificity benchmark.','Known-fusion recovery prior stays enabled except one planned comparison; do not extrapolate discovery sensitivity from BCR::ABL1.','Stale nonmaterial manifest prose says three full-reference cases; actual v1.1 list and runner have four, as reviewed here.','No clinical sensitivity percentage, limit of detection or patient fusion-negative conclusion is supported.']}
(P/'independent-design-review.json').write_text(json.dumps(review,indent=2)+'\n')
(P/'independent-plan-review.txt').write_text('APPROVED BEFORE EXECUTION\nManifest SHA256 '+manifest_hash+'\n\nIndependent reference, truth-pool, parent-control, full FASTQ and trimming audits all passed. Twenty paired datasets and 45 planned runs are approved for the bounded public-data question.\n\n'+'\n'.join('- '+x for x in review['limitations'])+'\n')
print(json.dumps({k:review[k] for k in ['status','manifest_sha256','source_pairs','eligible_pairs','mates_with_at_most_3_mismatches_checked','shortest_anchor_name']}))
