from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import csv,json,hashlib
O=Path('outputs/caris-followup/validation-package');records={r.id:str(r.seq) for r in SeqIO.parse(O/'edited-reference-cds.fasta','fasta')};vs=list(csv.DictReader((O/'variants.tsv').open(),delimiter='\t'));out=[]
for v in vs:
 g=v['gene'];rec=SeqIO.read(O/'references'/f'{v["refseq"]}.gb','genbank');feature=next(x for x in rec.features if x.type=='CDS');cds=str(feature.extract(rec.seq));assert cds==records[f'{g}_WT_edited_reference_CDS'];refprot=str(Seq(cds).translate(to_stop=True));assert refprot==feature.qualifiers['translation'][0]
 m=records[f'{g}_ALT_edited_reference_CDS'];mp=str(Seq(m[:len(m)//3*3]).translate(to_stop=True));assert len(mp)==int(v['conditional_mutant_aa_length'])
 if g=='BAP1':assert refprot[:56]==mp[:56] and mp[56:]=='KGLYLGG' and refprot[52:61]=='IEERRSRRK' and mp[52:61]=='IEERKGLYL'
 if g=='LATS1':assert mp==refprot[:126]
 if g=='LATS2':assert [i+1 for i,(a,b) in enumerate(zip(refprot,mp)) if a!=b]==[960] and refprot[959]=='G' and mp[959]=='R'
 out.append(dict(gene=g,refseq=v['refseq'],independent_Biopython_SeqIO_CDS_and_translation_match=True,mutant_aa_length=len(mp)))
for r in SeqIO.parse(O/'BAP1-peptide-nucleotide-haplotypes.fasta','fasta'):
 if '_ALT_' in r.id:
  expected='FKWIEERKGLYLGG*' if 'extended_context' in r.id else r.id.split('_')[1]
  assert str(r.seq.translate())==expected
out=dict(method='Separate Biopython GenBank feature extraction and translation compared to builder standalone codon table; full CDS protein must match NCBI /translation; peptide nucleotide FASTA translated independently',results=out)
(O/'independent-sequence-audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
