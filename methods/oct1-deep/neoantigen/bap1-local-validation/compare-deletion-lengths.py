"""Exact full peptide-window haplotype check of minor RNA 12bp deletion."""
from pathlib import Path
import json,collections
import pysam
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
SRC=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
prov=json.loads((OUT/'reference-provenance.json').read_text());win=prov['peptide_windows'][0]
st=win['start0'];en=win['end0'];wt=win['ref_nt']
templates={'reference':wt,'Caris_11bp_deletion':win['mutant_nt'],'minor_RNA_12bp_deletion':wt[:52408549-st]+wt[52408561-st:]}
out=[]
for kind in ['DNA','RNA']:
    with pysam.AlignmentFile(str(SRC/f'{kind}_TN26-279853.bam'),'rb',index_filename=str(ROOT/f'work/oct1-analysis/{kind}_TN26-279853.bam.bai')) as bam:
        rcount=collections.Counter();f=collections.defaultdict(set);c=collections.defaultdict(set);strands=collections.Counter()
        for r in bam.fetch('chr3',st,en):
            if r.flag&(4|256|512|1024|2048) or r.mapping_quality<20 or r.query_qualities is None:continue
            b={rp:q for q,rp in r.get_aligned_pairs(matches_only=True)}
            if st not in b or en-1 not in b:continue
            qa=b[st];qb=b[en-1]+1
            if min(r.query_qualities[qa:qb])<20:continue
            seq=r.query_sequence[qa:qb];call=next((k for k,v in templates.items() if seq==v),'other')
            rcount[call]+=1;k=(r.get_tag('RG') if r.has_tag('RG') else '',r.query_name);f[k].add(call);c[call].add((r.reference_start,r.reference_end,r.cigarstring,r.is_reverse,r.is_read1));strands[call+('_reverse' if r.is_reverse else '_forward')]+=1
        out.append({'kind':kind,'read_counts':dict(rcount),'fragment_counts':dict(collections.Counter(next(iter(v)) if len(v)==1 else 'discordant' for v in f.values())),'distinct_alignment_patterns':{k:len(v) for k,v in c.items()},'strands':dict(strands)})
(OUT/'minor-12bp-deletion-haplotype-check.json').write_text(json.dumps({'method':'Strict flags/MAPQ20, every nucleotide in the complete local nine-peptide coding window Q20; exact sequence matching; mate disagreement separated. Minor RNA deletion is a low-count observation, not a validated separate allele.','results':out},indent=2)+'\n')
print(json.dumps(out,indent=2))
