from pathlib import Path
from collections import Counter
import json,hashlib,pysam
P=Path(__file__).resolve().parent.parent
rc=lambda s:s.translate(str.maketrans('ACGTacgt','TGCAtgca'))[::-1]
fa=pysam.FastaFile('/Users/burhanazeem/.local/share/codex/caris-analysis/genome-fusion-reference/GRCh38.primary_assembly.genome.fa')
models=json.loads((P/'proposed-models.json').read_text());model=rc(fa.fetch('chr11',69651216,69651466))+rc(fa.fetch('chr3',125013072,125013572));assert all(x['sequence']==model for x in models)
shifts=[]
for d in range(-40,41):
 s=rc(fa.fetch('chr11',69651216-d,69651466))+rc(fa.fetch('chr3',125013072,125013572-d))
 assert len(s)==750
 if s==model:shifts.append(d)
marker=model[225:275];marker_rc=rc(marker)
pairs=json.loads((P/'pairs.json').read_text());result=[]
for p in pairs:
 levels={20:False,30:False};hits=[]
 for mate in [1,2]:
  seq=p[f'R{mate}_sequence'];qual=p[f'R{mate}_quality'];assert len(seq)==len(qual)
  for strand,pattern in [('+',marker),('-',marker_rc)]:
   start=0
   while True:
    loc=seq.find(pattern,start)
    if loc<0:break
    minq=min(ord(z)-33 for z in qual[loc:loc+50]);hits.append({'mate':mate,'orientation':strand,'start0':loc,'minQ':minq})
    for cutoff in levels:levels[cutoff]|=minq>=cutoff
    start=loc+1
 assert levels[20]==p['Q20_50nt_support']
 result.append({'pair_id':p['pair_id'],'family':p['sequence_family'],'Q20':levels[20],'Q30':levels[30],'fully_RC_mates':p['R1_sequence']==rc(p['R2_sequence']),'R1_length':len(p['R1_sequence']),'R2_length':len(p['R2_sequence']),'hits':hits})
q20=[p for p in result if p['Q20']]
out={'reference_models_independently_rebuilt':True,'two_proposed_models_identical':True,'model_sha256':hashlib.sha256(model.encode()).hexdigest(),'equivalent_junction_shifts_tested':[-40,40],'equivalent_junction_deltas':shifts,'delta_convention':'delta changes left RNA length250+delta and right RNA length500-delta; both genomic join coordinates decrease bydelta','Q20_pairs':len(q20),'Q30_pairs':sum(p['Q30'] for p in result),'Q20_fully_RC_pairs':sum(p['fully_RC_mates'] for p in q20),'Q20_sequence_families':len(set(p['family'] for p in q20)),'pair_details':result,'limits':'Exact overlapping mates/families are not independent molecules; reference-derived microhomology permits multiple breakpoint labels but does not distinguish biologicalRNA from library chimera.'}
assert shifts==list(range(-17,1));assert (out['Q20_pairs'],out['Q30_pairs'],out['Q20_fully_RC_pairs'],out['Q20_sequence_families'])==(11,9,10,7)
Path(__file__).with_name('independent-sequence-geometry.json').write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items() if k!='pair_details'})
