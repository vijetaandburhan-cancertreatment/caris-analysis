"""Exact-version RefSeq transcript consequences; predictions, not measured proteins."""
import hashlib, json, pathlib, time, csv
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
OUT=pathlib.Path(__file__).resolve().parent
SPEC=[
 {'gene':'BAP1','accession':'NM_004656.3','start_c':168,'end_c':178,'ref':'CCGGTCCCGGC','alt':'','reported_c':'c.168_178del11','reported_p':'p.R57fs','genomic_ref_deleted':'GCCGGGACCGG','genomic_to_transcript':'reverse_complement'},
 {'gene':'RASA1','accession':'NM_002890.2','start_c':747,'end_c':747,'ref':'G','alt':'','reported_c':'c.747delG','reported_p':'p.S250fs','genomic_ref_deleted':'G','genomic_to_transcript':'same_orientation'},
 {'gene':'APC','accession':'NM_000038.5','start_c':7730,'end_c':7730,'ref':'C','alt':'G','reported_c':'c.7730C>G','reported_p':'p.S2577*','genomic_ref_deleted':'C','genomic_to_transcript':'same_orientation'},
]
# Standard genetic code: independently specified table in T,C,A,G codon order.
AA='FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
BASE='TCAG'
TABLE={a+b+c:AA[i*16+j*4+k] for i,a in enumerate(BASE) for j,b in enumerate(BASE) for k,c in enumerate(BASE)}
NAMES={'R':'Arg','K':'Lys','S':'Ser','Q':'Gln','*':'Ter'}
def manual_translate(seq):
 ans=[]
 for i in range(0,len(seq)-2,3):
  amino=TABLE[seq[i:i+3]]
  if amino=='*':return ''.join(ans),i//3+1,seq[i:i+3]
  ans.append(amino)
 raise ValueError('No stop before transcript end')
def fasta(path,description,sequence):
 path.write_text('>'+description+'\n'+'\n'.join(sequence[i:i+70] for i in range(0,len(sequence),70))+'\n')
results=[]
GENOME_WINDOWS={x['gene']:x for x in json.loads((OUT.parent/'variant-read-support.hg38-reference-windows.json').read_text())}
for s in SPEC:
 gb=OUT/f'{s["accession"]}.gb';record=SeqIO.read(gb,'genbank')
 assert record.id==s['accession'],record.id
 features=[f for f in record.features if f.type=='CDS'];assert len(features)==1
 f=features[0];assert f.qualifiers.get('codon_start',['1'])==['1']
 assert f.location.strand==1 and len(f.location.parts)==1
 start=int(f.location.start);end=int(f.location.end);cds=str(f.extract(record.seq)).upper()
 assert len(cds)%3==0 and cds[:3]=='ATG'
 ref_protein=str(Seq(cds).translate(cds=True))
 assert ref_protein==f.qualifiers['translation'][0], 'reference CDS translation != NCBI protein'
 assert manual_translate(cds)[0]==ref_protein, 'independent normal translation mismatch'
 assert cds[s['start_c']-1:s['end_c']]==s['ref']
 gen_ref=s['genomic_ref_deleted']
 if s['genomic_to_transcript']=='reverse_complement':gen_ref=str(Seq(gen_ref).reverse_complement())
 assert gen_ref==s['ref'],'genomic/CDS allele mismatch'
 genomic_context_check=None
 if s['gene'] in GENOME_WINDOWS:
  w=GENOME_WINDOWS[s['gene']];gstart={'BAP1':52408550,'RASA1':87332560}[s['gene']]
  seq=w['dna'][gstart-20-w['start']:gstart+len(s['ref'])+20-w['start']].upper()
  if s['genomic_to_transcript']=='reverse_complement':seq=str(Seq(seq).reverse_complement())
  target=cds[s['start_c']-1-20:s['end_c']+20]
  assert seq==target,'extended genomic/transcript reference context mismatch'
  genomic_context_check={'public_hg38_context_matches_exact_RefSeq_CDS':True,'context_nt':len(seq),'genomic_deletion_start_0based':gstart,'flanks_each_side_nt':20,'source_url':w['source_url']}
 # Method A: edit extracted reference CDS, then translate using Biopython.
 edited=cds[:s['start_c']-1]+s['alt']+cds[s['end_c']:]
 a_protein=str(Seq(edited[:len(edited)//3*3]).translate(to_stop=True))
 # Method B: independently edit full mRNA at CDS-start+one-based c. offsets,
 # then translate manually from the original initiation site through first stop.
 tx=str(record.seq).upper();ts=start+s['start_c']-1;te=start+s['end_c']
 altered_tx=tx[:ts]+s['alt']+tx[te:]
 b_protein,stop_position,stop_codon=manual_translate(altered_tx[start:])
 assert a_protein==b_protein, 'independent mutant translation mismatch'
 assert stop_position<=len(edited)//3,'first stop not covered by edited reference CDS'
 first=next((i for i,(a,b) in enumerate(zip(ref_protein,a_protein)) if a!=b),min(len(ref_protein),len(a_protein)))
 prefix=ref_protein[:first];tail=a_protein[first:];assert a_protein[:first]==prefix
 newaa=a_protein[first] if first<len(a_protein) else '*'
 if len(s['alt'])!=len(s['ref']):
  fsTer=stop_position-(first+1)+1
  nomenclature=f'p.({NAMES[ref_protein[first]]}{first+1}{NAMES[newaa]}fsTer{fsTer})'
 else:
  fsTer=None;nomenclature=f'p.({NAMES[ref_protein[first]]}{first+1}{NAMES[newaa]})'
 orig_codon=cds[first*3:first*3+3];mut_codon=altered_tx[start+first*3:start+first*3+3]
 predicted_cds=altered_tx[start:start+stop_position*3]
 assert str(Seq(predicted_cds).translate(cds=True))==a_protein
 stem=f'{s["gene"]}_{s["accession"]}'
 fasta(OUT/f'{stem}.reference-cds.fasta',f'{s["accession"]} reference_CDS includes_terminal_stop',cds)
 fasta(OUT/f'{stem}.edited-reference-cds-span.fasta',f'{s["accession"]} {s["reported_c"]} edited_reference_CDS_span_not_reannotated_CDS',edited)
 fasta(OUT/f'{stem}.predicted-mutant-cds.fasta',f'{s["accession"]} {s["reported_c"]} PREDICTED_CDS_to_first_stop_including_stop',predicted_cds)
 fasta(OUT/f'{stem}.reference-protein.fasta',f'{f.qualifiers["protein_id"][0]} reference_protein',ref_protein)
 fasta(OUT/f'{stem}.predicted-mutant-protein.fasta',f'{s["accession"]} {s["reported_c"]} {nomenclature} PREDICTED_not_measured',a_protein)
 item={**s,'record_title':record.description,'record_date':record.annotations.get('date'),'record_sha256':hashlib.sha256(gb.read_bytes()).hexdigest(),'refseq_url':f'https://www.ncbi.nlm.nih.gov/nuccore/{record.id}','protein_accession':f.qualifiers['protein_id'][0],'cds_start_mrna_1based':start+1,'cds_end_mrna_1based_inclusive':end,'normal_mrna_nt':len(record.seq),'normal_cds_nt_including_stop':len(cds),'normal_protein_aa':len(ref_protein),'reference_deleted_or_replaced_nt':s['ref'],'first_changed_residue_1based':first+1,'reference_amino_acid':ref_protein[first],'predicted_new_amino_acid':newaa,'reference_codon_at_first_change':orig_codon,'mutant_codon_at_first_change':mut_codon,'predicted_protein_hgvs':nomenclature,'predicted_stop_codon_position_1based':stop_position,'predicted_stop_codon':stop_codon,'frameshift_Ter_distance':fsTer,'predicted_protein_aa':len(a_protein),'predicted_mutant_cds_nt_including_stop':len(predicted_cds),'unchanged_prefix_aa':first,'altered_tail_aa':len(tail),'altered_tail_sequence':tail,'reference_context_10aa_upstream':ref_protein[max(0,first-10):first+20],'mutant_context_10aa_upstream':a_protein[max(0,first-10):]+'*','independent_checks':{'exact_accession_and_version':True,'reference_cds_translation_matches_NCBI_qualifier':True,'manual_reference_translation_matches':True,'reference_allele_matches_c_coordinate':True,'genomic_reference_allele_matches_transcript_orientation':True,'independent_full_mRNA_edit_manual_translation_matches_CDS_edit_Biopython':True,'mutant_CDS_with_terminal_stop_translates_consistently':True}}
 item['extended_genomic_context_check']=genomic_context_check
 results.append(item)
report={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'biopython_version':Bio.__version__,'status':'computational predictions on named RefSeq isoforms; no demonstrated protein translation or antigen presentation','limits':['Local DNA and RNA support does not establish full-length transcript structure or phasing.','Other cis variants, alternative splicing, alternate initiation, RNA decay and protein turnover can change the outcome.','No matched-normal specimen establishes these alterations as somatic.','No HLA binding, peptide presentation, T-cell specificity or clinical efficacy was assessed.','Counts of altered residues are not counts of usable neoantigens.'], 'results':results}
(OUT/'transcript-consequences.json').write_text(json.dumps(report,indent=2)+'\n')
fields=['gene','accession','reported_c','reported_p','predicted_protein_hgvs','normal_cds_nt_including_stop','normal_protein_aa','predicted_protein_aa','first_changed_residue_1based','predicted_stop_codon_position_1based','altered_tail_aa','altered_tail_sequence','reference_context_10aa_upstream','mutant_context_10aa_upstream']
with (OUT/'transcript-consequences.tsv').open('w') as out:
 w=csv.DictWriter(out,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows({k:r[k] for k in fields} for r in results)
print(json.dumps(report,indent=2))
