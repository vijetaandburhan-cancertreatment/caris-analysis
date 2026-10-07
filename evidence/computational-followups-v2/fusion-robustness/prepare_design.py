"""Freeze a patient-independent public-fragment robustness benchmark before execution."""
from pathlib import Path
import gzip, json, hashlib, re, datetime, csv
import numpy as np
import pysam

P=Path(__file__).resolve().parent
C=Path('/Users/burhanazeem/.local/share/codex/caris-analysis')
PRE=C/'oct4-followup/genome-fusion'
REF=C/'genome-fusion-reference'
TOOLS=C/'genome-fusion-tools'
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def rc(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def fq(p):
    with gzip.open(p,'rt') as f:
        while h:=f.readline().rstrip():
            s=f.readline().rstrip();plus=f.readline().rstrip();q=f.readline().rstrip()
            assert h.startswith('@') and plus.startswith('+') and len(s)==len(q)
            yield {'header':h,'name':h[1:].split()[0].removesuffix('/1').removesuffix('/2'),'seq':s,'plus':plus,'qual':q}

fa=pysam.FastaFile(str(REF/'GRCh38.primary_assembly.genome.fa'))
spec={'BCR':{'tx':'ENST00000305877.13','chrom':'chr22','breakpoint1':23290413},'ABL1':{'tx':'ENST00000372348.7','chrom':'chr9','breakpoint1':130854064}}
for gene,d in spec.items():
    ex=[]
    for line in (REF/'compact_BCR_ABL1/annotation.gtf').read_text().splitlines():
        f=line.split('\t')
        if len(f)==9 and f[2]=='exon' and f'transcript_id "{d["tx"]}"' in f[8]:
            assert f[6]=='+';ex.append((int(f[3])-1,int(f[4])))
    ex=sorted(set(ex));assert ex
    d['exons0']=ex;d['sequence']=''.join(fa.fetch(d['chrom'],a,b).upper() for a,b in ex)
    # Breakpoints are the final retained BCR base and first retained ABL1 base.
    d['cut0']=sum(max(0,min(b,d['breakpoint1'] if gene=='BCR' else d['breakpoint1']-1)-a) for a,b in ex if a<(d['breakpoint1'] if gene=='BCR' else d['breakpoint1']-1))
fusion=spec['BCR']['sequence'][:spec['BCR']['cut0']]+spec['ABL1']['sequence'][spec['ABL1']['cut0']:]
j=spec['BCR']['cut0'];assert j>150 and len(fusion)-j>150
r1=list(fq(PRE/'public-controls/arriba.R1.fastq.gz'));r2=list(fq(PRE/'public-controls/arriba.R2.fastq.gz'));assert len(r1)==len(r2)==40
pairs={a['name']:(a,b) for a,b in zip(r1,r2)};assert len(pairs)==40 and all(a['name']==b['name'] for a,b in pairs.values())
assert len({(a['seq'],b['seq']) for a,b in pairs.values()})==40
truth=[];pool=[]
farray=np.frombuffer(fusion.encode(),dtype=np.uint8)
for name,(a,b) in pairs.items():
    mates=[]
    for mate,r in enumerate((a,b),1):
        n=len(r['seq']);windows=np.lib.stride_tricks.sliding_window_view(farray,n)
        trials=[]
        for orientation,s in [('+',r['seq']),('-',rc(r['seq']))]:
            distances=np.count_nonzero(windows!=np.frombuffer(s.encode(),dtype=np.uint8),axis=1)
            best=int(distances.min());starts=np.flatnonzero(distances==best).tolist()
            trials.append((best,orientation,starts,distances))
        best=min(t[0] for t in trials);loc=[(orient,start) for score,orient,starts,_ in trials if score==best for start in starts]
        unique=len(loc)==1;orient,start=loc[0]
        crossing=start<j<start+n
        left=j-start if crossing else None;right=start+n-j if crossing else None
        original_cut=(left if orient=='+' else right) if crossing else None
        # No gapped rescue or caller-defined support selection; ambiguous/minimally matching reads are excluded from pool.
        valid=unique and best<=3 and crossing and min(left,right)>=10
        mates.append({'mate':mate,'orientation':orient,'best_mismatches':best,'best_locations':len(loc),'start_on_fusion0':start,'end_on_fusion0':start+n,'crosses_junction':crossing,'left_anchor':left,'right_anchor':right,'original_read_junction_cut0':original_cut,'eligible_crossing_mate':valid})
    eligible=any(m['eligible_crossing_mate'] for m in mates)
    row={'name':name,'mates':mates,'eligible':eligible}
    if eligible:pool.append(name)
    truth.append(row)
assert len(pool)>=6,('Unexpected pool size',pool)

# A fixed TRUE parental background, not the non-spanning portion of a fusion library.
background={};parent_info=[]
source_names=list(pairs)
for gidx,(gene,d) in enumerate(spec.items()):
    seq=d['sequence'];cut=d['cut0']
    for k in range(16):
        length=[250,300,350,400][k%4]
        start=max(0,min(len(seq)-length,cut-450+k*47))
        name=f'UNFUSED_{gene}_{k:02d}'
        quality_source=source_names[gidx*16+k]
        qa,qb=pairs[quality_source]
        s1=seq[start:start+150];s2=rc(seq[start+length-150:start+length])
        assert len(s1)==len(s2)==150
        background[name]=({'header':'@'+name+'/1','name':name,'seq':s1,'plus':'+','qual':qa['qual']},{'header':'@'+name+'/2','name':name,'seq':s2,'plus':'+','qual':qb['qual']})
        parent_info.append({'name':name,'gene':gene,'transcript':d['tx'],'fragment_start0':start,'fragment_end0':start+length,'fragment_length':length,'quality_source_name':quality_source,'sequence_origin':'Exact normal GENCODE37 spliced parental transcript; sequence is newly synthetic; source qualities reused without modification.'})
assert len({(a['seq'],b['seq']) for a,b in background.values()})==32

datasets=[]
def write_dataset(label,selected,records,kind,extra=None):
    names=list(records);assert len(names)==len(set(names))
    seqpairs=[(a['seq'],b['seq']) for a,b in records.values()]
    assert len(seqpairs)==len(set(seqpairs)),'No duplicated identical fragment sequence pairs allowed'
    files=[]
    for m in range(2):
        path=P/'inputs'/f'{label}.R{m+1}.fastq.gz';path.parent.mkdir(exist_ok=True)
        raw=''.join('\n'.join([rs[m]['header'],rs[m]['seq'],rs[m]['plus'],rs[m]['qual']])+'\n' for rs in records.values()).encode()
        with path.open('wb') as out:
            with gzip.GzipFile(filename='',mode='wb',fileobj=out,mtime=0) as z:z.write(raw)
        files.append({'path':str(path),'sha256':sha(path),'uncompressed_sha256':hashlib.sha256(raw).hexdigest(),'bytes':path.stat().st_size})
    datasets.append({'id':label,'kind':kind,'source_fusion_names':selected,'input_pairs':len(names),'names':names,'files':files,**(extra or {})})

write_dataset('parental_zero',[],background,'synthetic unfused-parent zero-fusion control')
write_dataset('fixture_all40',list(pairs),pairs,'unaltered full public fixture; not independent of earlier control runs')
levels=[1,2,3,5]
for n in levels:
    seen=set()
    for seed in range(100):
        ordered=sorted(pool,key=lambda x:hashlib.sha256(f'fusion-robustness-v1:seed{seed}:{x}'.encode()).hexdigest())
        selected=sorted(ordered[:n]);key=tuple(selected)
        if key in seen:continue
        seen.add(key)
        label=f'support_{n:02d}_selection_{len(seen)}'
        records=dict(background);records.update({name:pairs[name] for name in selected})
        write_dataset(label,selected,records,'unmodified source pairs plus fixed synthetic parental background',{'seed':seed,'source_support_level':n})
        if len(seen)==3:break
    assert len(seen)==3
selected=sorted(pool);records=dict(background);records.update({n:pairs[n] for n in selected})
write_dataset('support_all',selected,records,'all eligible source crossing pairs plus fixed parental background',{'source_support_level':len(selected)})

# Prespecified short-anchor case: sequence-derived smallest eligible minor anchor (ties query name).
tr={r['name']:r for r in truth}
minor=lambda n:min(min(m['left_anchor'],m['right_anchor']) for m in tr[n]['mates'] if m['eligible_crossing_mate'])
short=sorted(pool,key=lambda n:(minor(n),n))[0]
records=dict(background);records[short]=pairs[short]
write_dataset('shortest_anchor', [short],records,'unaltered shortest-anchor source pair plus parental background',{'minor_anchor':minor(short)})

# Strong-source group chosen solely by sequence-derived crossing anchors; trim both sequence and qualities.
strong=sorted([n for n in pool if minor(n)>=40],key=lambda n:(-minor(n),n))[:4]
assert len(strong)==4,('Not four strong anchor source pairs',strong)
for cap in [None,12,20,30]:
    records=dict(background);changes=[]
    for name in strong:
        edited=[]
        for r,m in zip(pairs[name],tr[name]['mates']):
            a,b=0,len(r['seq'])
            if cap is not None and m['crosses_junction'] and m['best_locations']==1 and m['best_mismatches']<=3:
                cut=m['original_read_junction_cut0'];assert cut is not None
                if cut<=len(r['seq'])-cut:a=cut-cap
                else:b=cut+cap
                assert 0<=a<b<=len(r['seq'])
            new=dict(r);new['seq']=r['seq'][a:b];new['qual']=r['qual'][a:b]
            assert len(new['seq'])==len(new['qual'])
            edited.append(new)
            changes.append({'name':name,'mate':m['mate'],'original_start0_retained':a,'original_end0_retained':b,'retained_length':b-a,'original_junction_cut0':m['original_read_junction_cut0'],'new_junction_cut0':None if not m['crosses_junction'] else m['original_read_junction_cut0']-a,'sequence_and_quality_trimmed_identically':True})
        records[name]=tuple(edited)
    label='strong_original' if cap is None else f'strong_trim_anchor{cap}'
    write_dataset(label,strong,records,'unmodified strong anchors' if cap is None else 'end-trimmed public fragment pairs; no manufactured bases',{'anchor_cap':cap,'trim_records':changes})

files={'reference_fasta':REF/'GRCh38.primary_assembly.genome.fa','compact_gtf':REF/'compact_BCR_ABL1/annotation.gtf','STAR':TOOLS/'STAR_2.7.11b/MacOSX_x86_64/STAR','Arriba':TOOLS/'arriba_v2.5.1/arriba','source_R1':PRE/'public-controls/arriba.R1.fastq.gz','source_R2':PRE/'public-controls/arriba.R2.fastq.gz'}
# Do not scan the 3 GB unchanged reference again; exact manifest already verified in prior run.
provenance={k:{'path':str(p),'sha256':sha(p) if k!='reference_fasta' else None,'bytes':p.stat().st_size} for k,p in files.items()}
provenance['reference_fasta']['hash_source']='Existing official reference and prior completed pipeline manifest; unchanged resident FASTA.'
out={'version':'v1.1','frozen_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Public-fixture pipeline robustness, not clinical sensitivity, limit of detection, biological reproducibility, patient validation or generalizable percentage.','truth_method':'Build fixed spliced BCR::ABL1 from declared GENCODE37 transcript breakpoints; ungapped Hamming alignment of each original read and reverse complement; eligible pair has at least one uniquely best mate with <=3 mismatches and >=10 nt on each junction side. Selection does not use new or previous STAR/Arriba acceptance. All source read classification retained, including exclusions.','junction':{'BCR_last_retained_base1':'chr22:23290413','ABL1_first_retained_base1':'chr9:130854064','fusion_transcript_boundary0':j,'fusion_transcript_length':len(fusion),'fusion_transcript_sha256':hashlib.sha256(fusion.encode()).hexdigest()},'parent_transcripts':{g:{k:v for k,v in d.items() if k!='sequence'}|{'sequence_sha256':hashlib.sha256(d['sequence'].encode()).hexdigest()} for g,d in spec.items()},'source_provenance':provenance,'truth_reads':truth,'eligible_source_names':pool,'shortest_anchor_name':short,'strong_anchor_names':strong,'parent_fragments':parent_info,'datasets':datasets,'compact_matrix':[{'dataset':d['id'],'sparseD':s} for d in datasets for s in [1,8]],'full_reference_D8_transfer_checks':['parental_zero','support_all','strong_original','strong_trim_anchor12'],'no_prior_comparison':{'dataset':'support_01_selection_1','sparseD':1,'omit_arriba_flags':['-k'],'keep_tags_flag_t':True,'reason':'One predeclared check of known-fusion recovery prior; no filter thresholds changed.'},'resources':{'max_process_RSS_GiB':12,'disk_min_GiB':3,'STAR_threads':2,'execution':'Sequential STAR to small BAM, then Arriba; gzip all caller logs during execution to avoid repetitive unknown-contig warning disk usage; no new index/download.','expected_max_runs':2*len(datasets)+5},'outcomes':['Exact breakpoints accepted/discarded, confidence and filter reasons','Distinct input query names vs distinct STAR target-junction names vs Arriba supporting names','Per-name target/competing-locus alignment geometry, CIGAR, anchors, MAPQ, NH and NM','Target support identity intersections/differences D1 versus D8','Warnings and resource metrics; no detection-rate extrapolation'],'fixed_settings':'STAR2.7.11b, Arriba2.5.1, GENCODE37; identical declared mapping/filter parameters; SAindex12 and overhang160; only sparseD varies in compact pairwise comparisons. Parent background fixed across support/anchor datasets. Source pair identity/sequences/quality unchanged except explicitly trimmed datasets. No duplicate copies of identical source fragments.','limitations':['Selected subsets overlap and share one synthetic fixture; not independent biological replicates.','Ungapped sequence-derived eligibility threshold is a benchmark definition, not all possible true fragment truth.','Compact reference lacks competing loci; full-reference D8 transfer checks address only three prespecified cases.','Known BCR::ABL1 prior stays enabled except one declared no-prior comparison.','Trimming also changes read length and mate geometry, so effects cannot be assigned solely to anchor length.','Synthetic parental controls reuse fixture quality strings and sample a tiny portion of normal transcript diversity.','No patient data, new clinical claims, whole-RNA rerun or changes to earlier outputs.']}
(P/'design-manifest.json').write_text(json.dumps(out,indent=2))
(P/'truth-fusion-transcript.fa').write_text('>GENCODE37_BCR_ABL1_fixed_junction\n'+fusion+'\n')
with (P/'normal-parent-transcripts.fa').open('w') as f:
    for g,d in spec.items():f.write('>'+g+'|'+d['tx']+'\n'+d['sequence']+'\n')
(P/'DESIGN.txt').write_text('Prespecified public fusion robustness benchmark v1\n\n'+out['purpose']+'\n\n'+out['truth_method']+'\n\n'+out['fixed_settings']+'\n\n'+f'{len(datasets)} datasets, each run on compact D1 and D8; four declared full-reference D8 transfer checks and one no-known-fusion-recovery comparison.\n\n'+'\n'.join('- '+x for x in out['limitations'])+'\n')
(P/'design-freeze.json').write_text(json.dumps({'manifest_sha256':sha(P/'design-manifest.json'),'script_sha256':sha(Path(__file__)),'frozen_utc':out['frozen_utc'],'status':'Awaiting independent design review; no benchmark alignment executed.'},indent=2))
print(json.dumps({'eligible_source_names':pool,'shortest':short,'strong':strong,'datasets':len(datasets),'planned_runs':out['resources']['expected_max_runs'],'freeze':json.loads((P/'design-freeze.json').read_text())},indent=2))
