#!/usr/bin/env python3
"""Bounded Hamming-distance reference context; not genome-wide approximate mapping."""
import json,pathlib,pysam
base=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis/oct4-fusion-read-audit');a=json.load(open(base/'patient/junction-read-audit.json'));path='/Users/burhanazeem/.local/share/codex/caris-analysis/genome-fusion-reference/GRCh38.primary_assembly.genome.fa';fa=pysam.FastaFile(path);out=[]
def rc(s):return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
for i in [2,20,23,25,26,27]:
 r=a['rows'][i];p=next(x['sequence'] for x in r['junction_patterns'] if x['arm_bases']==25);rec={'candidate_id':r['candidate_id'],'gene1':r['gene1'],'gene2':r['gene2'],'marker':p,'loci':[]}
 for side in [1,2]:
  ch,pos=r[f'breakpoint{side}'].split(':');pos=int(pos);reg=fa.fetch(ch,pos-150,pos+150).upper();hits=[]
  for ori,seq in [('forward',reg),('reverse',rc(reg))]:
   for start in range(len(seq)-len(p)+1):
    diff=[j for j,(x,y) in enumerate(zip(seq[start:start+len(p)],p)) if x!=y]
    if len(diff)<=2:hits.append({'orientation':ori,'oriented_window_start0':start,'mismatch_positions0':diff,'normal_sequence':seq[start:start+len(p)]})
  rec['loci'].append({'side':side,'contig':ch,'reference_window_start0':pos-150,'reference_window_end0':pos+150,'hits_with_at_most2_mismatches':hits})
 out.append(rec)
(base/'patient-evidence/normal-locus-neighbor-context.json').write_text(json.dumps({'reference':path,'rows':out,'limits':['Only300nt around each reported genomic locus, ungapped at most2 mismatches. Not a global approximate-alignment or variant-calling method.','Near-normal context supports an alternate explanation but cannot by itself prove an artifact or identify the actual allele.']},indent=2)+'\n')
