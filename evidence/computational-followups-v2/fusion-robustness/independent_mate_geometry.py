"""Independent physical-pair geometry of frozen eligible public source fragments."""
from pathlib import Path
import gzip,json,re,hashlib,datetime
P=Path(__file__).resolve().parent
def rc(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def records(p):
    with gzip.open(p,'rt') as f:
        while h:=f.readline():
            s=f.readline().strip();f.readline();q=f.readline().strip()
            yield re.sub(r'/[12]$','',h[1:].split()[0]),s,q
def bestplacements(query,ref):
    # Exhaustive seeded search for every ungapped <=3-mismatch alignment.
    found=[];n=len(query)
    for ori,seq in [('+',query),('-',rc(query))]:
        positions=set()
        for k in range(4):
            a,b=n*k//4,n*(k+1)//4;at=ref.find(seq[a:b])
            while at>=0:
                if 0<=at-a<=len(ref)-n:positions.add(at-a)
                at=ref.find(seq[a:b],at+1)
        for start in positions:
            mm=sum(a!=b for a,b in zip(seq,ref[start:start+n]))
            if mm<=3:found.append((mm,ori,start))
    if not found:return []
    low=min(x[0] for x in found);return sorted(x for x in found if x[0]==low)
d=json.loads((P/'design-manifest.json').read_text());ref=''.join((P/'truth-fusion-transcript.fa').read_text().splitlines()[1:]);j=d['junction']['fusion_transcript_boundary0']
src=[{n:(s,q) for n,s,q in records(d['source_provenance'][f'source_R{i}']['path'])} for i in [1,2]]
rows=[]
for name in d['eligible_source_names']:
    mates=[]
    for i in range(2):
        seq,q=src[i][name];hits=bestplacements(seq,ref)
        item={'mate':i+1,'best_placements_at_most_3_mismatches':hits,'length':len(seq),'uniquely_placed':len(hits)==1}
        if len(hits)==1:
            mm,ori,start=hits[0];aligned=seq if ori=='+' else rc(seq);qual=q if ori=='+' else q[::-1]
            mismatches=[{'construct_base0':start+k,'original_read_base0':k if ori=='+' else len(seq)-1-k,'base_quality':ord(qual[k])-33} for k,(x,y) in enumerate(zip(aligned,ref[start:start+len(seq)])) if x!=y]
            item.update(orientation=ori,construct_start0=start,construct_end0=start+len(seq),mismatches=mm,mismatch_positions=mismatches,crosses_junction=start<j<start+len(seq),left_anchor=j-start if start<j<start+len(seq) else None,right_anchor=start+len(seq)-j if start<j<start+len(seq) else None)
        mates.append(item)
    row={'name':name,'mates':mates,'strong_group':name in d['strong_anchor_names'],'shortest_case':name==d['shortest_anchor_name']}
    if all(x['uniquely_placed'] for x in mates):
        starts=[x['construct_start0'] for x in mates];ends=[x['construct_end0'] for x in mates];orient=[x['orientation'] for x in mates]
        opposite=orient[0]!=orient[1];plus=next((x for x in mates if x['orientation']=='+'),None);minus=next((x for x in mates if x['orientation']=='-'),None)
        inward=opposite and plus['construct_start0']<=minus['construct_start0'] and plus['construct_end0']<=minus['construct_end0']
        a,b=max(starts),min(ends);conflicts=[]
        if a<b:
            orienteds=[]
            for i,m in enumerate(mates):
                seq,q=src[i][name];orienteds.append((seq,q) if m['orientation']=='+' else (rc(seq),q[::-1]))
            for c in range(a,b):
                ia,ib=c-starts[0],c-starts[1]
                if orienteds[0][0][ia]!=orienteds[1][0][ib]:
                    conflicts.append({'construct_base0':c,'mate1_quality':ord(orienteds[0][1][ia])-33,'mate2_quality':ord(orienteds[1][1][ib])-33})
        row.update(opposite_orientations=opposite,inward_geometry_under_construct=inward,outer_fragment_start0=min(starts),outer_fragment_end0=max(ends),outer_fragment_span=max(ends)-min(starts),overlap_bases=max(0,b-a),overlap_disagreements=len(conflicts),overlap_disagreement_positions=conflicts,construct_based_overlap_disagreement_fraction=len(conflicts)/(b-a) if b>a else None)
    rows.append(row)
result={'status':'PASS','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'manifest_sha256':hashlib.sha256((P/'design-manifest.json').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'rows':rows,'interpretation':['Frozen eligibility required at least one uniquely placed crossing mate; this additional audit checks both mates separately, without altering selection.','Inwardness and plausible insert length on the fixed fusion construct support a possible physical pair, not molecular independence or genome-wide mappability.','An overlap mismatch may reflect an ordinary sequencing error; this audit cannot assign the cause of a STAR rejection.','No reads were removed or altered because of these post-freeze checks.']}
default_source=P.parent/'oct4-followup/genome-fusion/STAR-parametersDefault'
assert re.search(r'^peOverlapMMp\s+0.01\s*$',default_source.read_text(),re.M)
result['merging_hypothesis']={'source_path':str(default_source),'source_sha256':hashlib.sha256(default_source.read_bytes()).hexdigest(),'documented_default_peOverlapMMp':0.01,'actual_run_peOverlapNbasesMin':10,'construct_overlap_fraction_above_default':[r['name'] for r in rows if (r['construct_based_overlap_disagreement_fraction'] or 0)>0.01],'caveat':'These fractions are independently derived from a fixed construct alignment. Being above the documented maximum overlap mismatch proportion is a possible reason not to merge the mates, not a demonstrated causal explanation for STAR rejection. Pair28 is recovered by compact D1 as an unmerged chimera. No parameter was altered to rescue a result.'}
(P/'independent-mate-geometry.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'pairs':len(rows),'both_unique':sum(all(m['uniquely_placed'] for m in r['mates']) for r in rows),'opposed_inward':sum(r.get('inward_geometry_under_construct',False) for r in rows),'span_range':[min(r['outer_fragment_span'] for r in rows),max(r['outer_fragment_span'] for r in rows)],'unmapped_strong_74':next(r for r in rows if r['name']=='BCR-ABL1-74')},indent=2))
