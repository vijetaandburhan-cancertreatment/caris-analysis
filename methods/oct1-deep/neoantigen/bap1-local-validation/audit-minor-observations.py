"""Quality audit of low-fraction nearby observations; no new variant calls."""
from pathlib import Path
import collections,csv,json,statistics
import pysam

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
SRC=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
OLD=ROOT/'work/oct1-analysis'
x=json.loads((OUT/'results.json').read_text())
snvs=[z for r in x['results'] for z in r['nearby_SNV_candidates_above_threshold']]
indels=list(csv.DictReader((OUT/'local-indel-observations.tsv').open(),delimiter='\t'))
results=[]
for kind in ['DNA','RNA']:
    with pysam.AlignmentFile(str(SRC/f'{kind}_TN26-279853.bam'),'rb',index_filename=str(OLD/f'{kind}_TN26-279853.bam.bai')) as bam:
        candidates=[{'type':'SNV','p0':z['pos1']-1,'alt':z['alt_forward'],'coding_pos1':z['coding_pos1']} for z in snvs if z['kind']==kind]
        candidates += [{'type':'indel','p0':int(z['position0']),'op':2 if z['operation']=='D' else 1,'length':int(z['length']),'coding_pos1':z['coding_pos1_at_start']} for z in indels if z['kind']==kind and z['is_equivalent_known_BAP1_deletion']=='False']
        for v in candidates:
            p=v['p0'];detail=[];fragments=set();clean=set();coords=set();local_cigars=collections.Counter()
            for r in bam.fetch('chr3',p,p+1):
                if r.flag&(4|256|512|1024|2048) or r.mapping_quality<20 or r.query_qualities is None:continue
                b={rp:q for q,rp in r.get_aligned_pairs(matches_only=True)}
                if v['type']=='SNV':
                    q=b.get(p)
                    if q is None or r.query_qualities[q]<20 or r.query_sequence[q]!=v['alt']:continue
                else:
                    rp=r.reference_start;qp=0;found=False
                    for op,n in r.cigartuples:
                        if op==v['op'] and rp==p and n==v['length']:q=qp;found=True;break
                        if op in (0,7,8):rp+=n;qp+=n
                        elif op in (2,3):rp+=n
                        elif op in (1,4):qp+=n
                    if not found:continue
                    pp=list(range(p-3,p))+list(range(p+(v['length'] if v['op']==2 else 0),p+(v['length'] if v['op']==2 else 0)+3))
                    if not all(z in b and r.query_qualities[b[z]]>=20 for z in pp):continue
                k=(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name);fragments.add(k)
                coords.add((r.reference_start,r.reference_end,r.cigarstring,r.is_reverse,r.is_read1))
                soft=any(op==4 for op,n in r.cigartuples);distance=min(q-r.query_alignment_start,r.query_alignment_end-1-q)
                if not soft and distance>=5:clean.add(k)
                detail.append({'reverse':r.is_reverse,'softclipped':soft,'distance_to_alignment_end':distance,'BQ':int(r.query_qualities[min(q,len(r.query_qualities)-1)]),'MAPQ':r.mapping_quality})
                # Read sequence and names are deliberately not written.
                local_cigars[r.cigarstring]+=1
            results.append({'kind':kind,**v,'reads':len(detail),'query_name_fragments':len(fragments),'distinct_alignment_patterns':len(coords),'forward_reads':sum(not z['reverse'] for z in detail),'reverse_reads':sum(z['reverse'] for z in detail),'softclipped_reads':sum(z['softclipped'] for z in detail),'near_end_less_5_reads':sum(z['distance_to_alignment_end']<5 for z in detail),'clean_fragments':len(clean),'median_end_distance':statistics.median(z['distance_to_alignment_end'] for z in detail) if detail else None,'median_base_quality':statistics.median(z['BQ'] for z in detail) if detail else None,'CIGAR_counts':dict(local_cigars)})
(OUT/'minor-observation-audit.json').write_text(json.dumps({'scope':'Low-fraction observations retained for artifact review, not validated variant calls. RNA fragment support may be PCR-duplicated; coordinate diversity recorded.','results':results},indent=2)+'\n')
print(json.dumps(results,indent=2))
