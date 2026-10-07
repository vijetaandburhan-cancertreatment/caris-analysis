"""Reference-based conditional consequence only, not protein validation."""
import json,pathlib
from Bio.Seq import Seq
B=pathlib.Path(__file__).resolve().parent
ann=json.loads((B/'target-annotation.json').read_text())['genes']
x=json.loads((B/'lats-exploratory-screen.json').read_text())
for v in x['candidates']:
 g=ann[v['gene']];r=json.loads((B/f'{v["gene"]}.hg38-reference.json').read_text());seq=r['dna'].upper();start=g['gene_start0']
 genomic_cds=''.join(seq[s-start:e-start] for s,e in g['selected_CDS']);p=v['pos1']-1;off=0;inside=False
 for s,e in g['selected_CDS']:
  if s<=p<e:off+=p-s;inside=True;break
  off+=e-s
 if not inside:v['conditional_protein']='flank_not_CDS';continue
 assert genomic_cds[off:off+len(v['ref'])]==v['ref']
 mutant_genomic_cds=genomic_cds[:off]+v['alt']+genomic_cds[off+len(v['ref']):]
 cds=genomic_cds if g['strand']=='+' else str(Seq(genomic_cds).reverse_complement())
 alt=mutant_genomic_cds if g['strand']=='+' else str(Seq(mutant_genomic_cds).reverse_complement())
 rp=str(Seq(cds).translate());ap=str(Seq(alt[:len(alt)//3*3]).translate());assert rp.startswith('M') and '*' not in rp
 diff=next((i for i,(a,b) in enumerate(zip(rp,ap)) if a!=b),None)
 v.update(transcript=g['selected_transcript'],reference_CDS_length=len(cds),reference_protein_length=len(rp),conditional_protein=('synonymous' if diff is None else rp[diff]+str(diff+1)+ap[diff]+(' (1bp insertion; immediate stop)' if len(cds)!=len(alt) and ap[diff]=='*' else 'fs' if len(cds)!=len(alt) else '')),first_altered_aa_1based=diff+1 if diff is not None else None)
 if len(v['ref'])==len(v['alt'])==1:
  cp=off if g['strand']=='+' else len(cds)-off-1
  refbase=cds[cp];altbase=alt[cp];v['conditional_coding']=f'c.{cp+1}{refbase}>{altbase}'
 else:
  # Find all equivalent 1base transcript insertions, choose HGVS 3-prime most.
  assert len(alt)==len(cds)+1
  possible=[(k,alt[k]) for k in range(len(alt)) if cds[:k]+alt[k]+cds[k:]==alt]
  k,b=max(possible);v['transcript_insertion_offset0']=k;v['transcript_inserted_base']=b
  v['conditional_coding']=f'c.{k}dup' if k>0 and cds[k-1]==b else f'c.{k}_{k+1}ins{b}'
  v['equivalent_transcript_insertion_offsets0']=[k for k,b in possible]
  v['first_stop_aa_1based']=ap.index('*')+1 if '*' in ap else None
  v['novel_amino_acid_tail_before_first_stop']=0 if diff is not None and ap[diff]=='*' else None
  run=seq[p-start+1];lo=p-start+1;hi=lo+1
  while lo>0 and seq[lo-1]==run:lo-=1
  while hi<len(seq) and seq[hi]==run:hi+=1
  v['insertion_adjacent_reference_homopolymer']={'base':run,'length':hi-lo,'start0':start+lo,'end0':start+hi}
  # Reference-insertion representations equivalence near primary CIGAR anchor.
  wst=p-40;wend=p+41;wr=seq[wst-start:wend-start];wa=wr[:p-wst]+v['alt']+wr[p-wst+len(v['ref']):]
  positions=[wst+k for k in range(len(wr)+1) if wr[:k]+v['alt'][len(v['ref']):]+wr[k:]==wa]
  v['equivalent_genomic_insertion_boundaries0']=positions
x['annotation_limit']='GENCODE v37 MANE CDS assembled from public hg38 reference; exact named transcript only. RNA read support does not verify full transcript/protein. No somatic status, pathogenic classification, or clinical actionability established.'
(B/'lats-conditional-annotation.json').write_text(json.dumps(x,indent=2)+'\n')
for v in x['candidates']:print(v['gene'],v['pos1'],v['conditional_coding'],v['conditional_protein'],v.get('equivalent_genomic_insertion_boundaries0'))
