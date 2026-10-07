"""Allele-specific consequence of the nearby NLRC3 SNV; no remote submission."""
import json,gzip
from pathlib import Path
from Bio.Seq import Seq
P=Path(__file__).resolve().parent;ROOT=P.parents[2];tx='ENST00000359128.10'
a=next(r for r in json.loads((P/'selected-coding-transcripts.json').read_text()) if r['transcript']==tx)
r=next(r for r in json.loads((P/'conditional-transcript-consequences.json').read_text())['results'] if r['gene']=='NLRC3')
positions=[p for s,e in sorted(a['CDS'],reverse=True) for p in range(e-1,s-1,-1)];i=positions.index(3563596);ex=[]
with gzip.open(ROOT/'work/oct1-deep/genomics/gencode.v37.annotation.gtf.gz','rt') as f:
 for line in f:
  if line.startswith('#'):continue
  a=line.rstrip().split('\t')
  if a[2]=='exon' and f'transcript_id "{tx}"' in a[8]:ex.append((int(a[3])-1,int(a[4])))
mrna=[p for s,e in sorted(ex,reverse=True) for p in range(e-1,s-1,-1)];offset=mrna.index(positions[0]);seq=json.loads((P/'selected-reference-transcripts.json').read_text())[tx]['sequence'];coding=seq[offset:offset+len(positions)];assert str(Seq(coding).translate())==r['reference_protein'];s=list(coding);assert s[i]=='C';s[i]='T';new=''.join(s);newprotein=str(Seq(new).translate(to_stop=True))
obj={'gene':'NLRC3','transcript':tx,'variant':'chr16:3563597G>A','coding_change':f'c.{i+1}C>T','protein_change':f'p.{r["reference_protein"][i//3]}{i//3+1}{newprotein[i//3]}','reference_codon':coding[3*(i//3):3*(i//3)+3],'alternate_codon':new[3*(i//3):3*(i//3)+3],'conditional_mutant_length':len(newprotein),'public_overlap_id':'rs773155391','note':'Nearby SNV is observed on reads without the insertion. Local read phasing only: not a matched normal comparator, biallelic tumor inference, pathogenicity result or clinical actionability claim. Public record-wide worst consequence stop_gained must not be applied to this observed G>A allele; selected transcript consequence is Ser447Leu.'}
(P/'audit-nlrc3-neighbor-consequence.json').write_text(json.dumps(obj,indent=2)+'\n');print(obj)
