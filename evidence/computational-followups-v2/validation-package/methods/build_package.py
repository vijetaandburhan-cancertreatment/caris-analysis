#!/usr/bin/env python3
"""Build a local, research-only BAP1/LATS sequence handoff from pinned evidence.
No sequence is transmitted; only two named public RefSeq records were retrieved.
"""
import csv,hashlib,json,re,shutil,textwrap
from pathlib import Path
ROOT=Path.cwd();W=ROOT/'work/oct4-followup/validation-package';O=ROOT/'outputs/caris-followup/validation-package';O.mkdir(parents=True,exist_ok=True)
C=str.maketrans('ACGT','TGCA')
COD=dict(zip([a+b+c for a in 'TCAG' for b in 'TCAG' for c in 'TCAG'],'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'))
rc=lambda s:s.translate(C)[::-1]
def tr(s):return ''.join(COD[s[i:i+3]] for i in range(0,len(s)-2,3))
def j(p):return json.loads((ROOT/p).read_text())
def tsv(name,rows):
 with (O/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def fasta(name,rows):
 with (O/name).open('w') as f:
  for r in rows:f.write('>'+r['id']+' '+r['description']+'\n'+'\n'.join(textwrap.wrap(r['sequence'],70))+'\n')
def gb(p):
 t=p.read_text();s=''.join(re.findall('[acgt]+',t.split('ORIGIN')[1].split('//')[0])).upper();a,b=map(int,re.search(r'\n     CDS\s+(\d+)\.\.(\d+)',t).groups());return s[a-1:b],re.search(r'/protein_id="([^"]+)"',t).group(1)
sources=set(); checks=[]
variants=[dict(gene='BAP1',chrom='chr3',accession='NC_000003.12',pos1=52408550,ref='CGCCGGGACCGG',alt='C',refseq='NM_004656.3',enst='ENST00000460680.6',coding='c.168_178del',protein='p.(Arg57LysfsTer8)',g_hgvs='g.52408551_52408561del',status='reported deletion; research peptide hypothesis'),dict(gene='LATS1',chrom='chr6',accession='NC_000006.12',pos1=149695191,ref='T',alt='TA',refseq='NM_004690.4',enst='ENST00000543571.6',coding='c.378dup',protein='p.(Asn127Ter)',g_hgvs='g.149695192dup',status='exploratory unreported insertion; orthogonal confirmation required'),dict(gene='LATS2',chrom='chr13',accession='NC_000013.11',pos1=20975259,ref='C',alt='T',refseq='NM_014572.3',enst='ENST00000382592.5',coding='c.2878G>A',protein='p.(Gly960Arg)',g_hgvs='g.20975259C>T',status='exploratory unreported missense; function unknown')]
ann=j('work/oct1-deep/genomics/target-annotation.json')['genes'];genrows=[];transrows=[];protrows=[];allcds=[]
for v in variants:
 gene=v['gene'];rp='work/oct1-deep/splicing/reference/BAP1.hg38-reference.json' if gene=='BAP1' else f'work/oct1-deep/genomics/{gene}.hg38-reference.json';sources.add(rp);r=j(rp);seq=r['dna'].upper();start=r['start'];p=v['pos1']-1;off=p-start;assert seq[off:off+len(v['ref'])]==v['ref']
 a=p-200;b=p+len(v['ref'])+200;wt=seq[a-start:b-start];alt=wt[:200]+v['alt']+wt[200+len(v['ref']):]
 for role,s in [('WT',wt),('ALT',alt)]:genrows.append(dict(id=f'{gene}_{role}_genomic_plus',description=f'GRCh38 {v["chrom"]}:{a+1}-{b}; plus strand; ALT coordinates anchored to WT; reference-derived single edit; not an oligo',sequence=s))
 v.update(build='GRCh38',strand='-',genomic_window_start1=a+1,genomic_window_end1=b)
 gbpath=ROOT/f'work/oct1-analysis/transcripts/{v["refseq"]}.gb' if gene=='BAP1' else W/'references'/f'{v["refseq"]}.gb';sources.add(str(gbpath.relative_to(ROOT)));cds,pr=gb(gbpath);v['protein_accession']=pr
 if gene=='BAP1':
  assert cds[167:178]=='CCGGTCCCGGC';m=cds[:167]+cds[178:];cs=124;ce=223;aa=57
  # The proven target exon matches the public genome and pins transcript orientation.
  prf=j('work/oct1-deep/neoantigen/bap1-local-validation/reference-provenance.json');ex=prf['target_exon'];assert rc(seq[ex['start0']-start:ex['end0']-start])==cds[ex['coding_start1']-1:ex['coding_end1']]
  tw=cds[cs-1:ce];tm=cds[cs-1:167]+cds[178:ce]
 elif gene=='LATS1':
  gcds=''.join(seq[x-start:y-start] for x,y in ann[gene]['selected_CDS']);assert rc(gcds)==cds[:-3]
  assert cds[377]=='T';m=cds[:378]+'T'+cds[378:];cs=334;ce=423;aa=127;tw=cds[cs-1:ce];tm=cds[cs-1:378]+'T'+cds[378:ce]
  # Separate genomic allele application must reproduce the same transcript edit.
  coords=[x for s,e in ann[gene]['selected_CDS'] for x in range(s,e)];at=coords.index(p);ga=gcds[:at]+v['alt']+gcds[at+len(v['ref']):];assert rc(ga)==m[:-3]
 else:
  gcds=''.join(seq[x-start:y-start] for x,y in ann[gene]['selected_CDS']);assert rc(gcds)==cds[:-3]
  assert cds[2877]=='G';m=cds[:2877]+'A'+cds[2878:];cs=2833;ce=2925;aa=960;tw=cds[cs-1:ce];tm=m[cs-1:ce]
  coords=[x for s,e in ann[gene]['selected_CDS'] for x in range(s,e)];at=coords.index(p);ga=gcds[:at]+v['alt']+gcds[at+len(v['ref']):];assert rc(ga)==m[:-3]
 wtpr=tr(cds).split('*')[0];mp=tr(m).split('*')[0];v.update(wt_aa_length=len(wtpr),conditional_mutant_aa_length=len(mp),cds_window_ref_start1=cs,cds_window_ref_end1=ce)
 if gene=='BAP1':assert len(mp)==63 and mp[56:]=='KGLYLGG'
 if gene=='LATS1':assert len(mp)==126 and mp==wtpr[:126]
 if gene=='LATS2':assert len(mp)==1088 and wtpr[959]=='G' and mp[959]=='R' and sum(x!=y for x,y in zip(wtpr,mp))==1
 for role,s in [('WT',tw),('ALT',tm)]:transrows.append(dict(id=f'{gene}_{role}_CDS_window',description=f'{v["refseq"]} transcript sense c.{cs}_{ce} reference span; DNA alphabet T; ALT length varies with edit; not patient full haplotype',sequence=s))
 lo=max(0,aa-26);hi=min(len(wtpr),aa+24)
 for role,s in [('WT',wtpr[lo:hi]),('ALT',mp[lo:min(hi,len(mp))])]:protrows.append(dict(id=f'{gene}_{role}_protein_window',description=f'{pr} aa {lo+1}-{min(hi,len(wtpr) if role=="WT" else len(mp))}; {"actual predicted termination follows immediately" if role=="ALT" and len(mp)<hi else "window only, protein continues"}',sequence=s))
 for role,s in [('WT',cds),('ALT',m)]:allcds.append(dict(id=f'{gene}_{role}_edited_reference_CDS',description=f'{v["refseq"]}; full reference CDS span edited; includes stop and may include untranslated sequence after premature stop; not synthesis construct',sequence=s))
 checks.append(dict(gene=gene,refseq=v['refseq'],genomic_ref_allele_exact=True,refseq_target_context_matches_GRCh38=True,conditional_consequence=v['protein'],wt_aa_length=len(wtpr),mutant_aa_length=len(mp)))
fasta('genomic-windows.fasta',genrows);fasta('transcript-windows.fasta',transrows);fasta('protein-windows.fasta',protrows);fasta('edited-reference-cds.fasta',allcds);tsv('variants.tsv',variants)
tsv('sequence-windows.tsv',[dict(molecule=m,**r) for m,rs in [('genomic_DNA_plus',genrows),('transcript_CDS_sense',transrows),('protein',protrows)] for r in rs])
# Exact peptide haplotypes (both strand conventions explicitly included).
pref='work/oct1-deep/neoantigen/bap1-local-validation/reference-provenance.json';sources.add(pref);peps=[]
for w in j(pref)['peptide_windows']:
 for role,key in [('WT','ref_nt'),('ALT','mutant_nt')]:
  s=rc(w[key]);peps.append(dict(id=f'BAP1_{w["label"]}_{role}_CDS_sense',description=f'GRCh38 minus-strand source chr3:{w["start0"]+1}-{w["end0"]}; reference span c.{w["ref_c_start"]}_{w["ref_c_end"]}; ALT is actual translated reading frame',sequence=s))
  if role=='ALT':assert tr(s)==w['expected_peptide']
fasta('BAP1-peptide-nucleotide-haplotypes.fasta',peps)
# Existing independent read support; do not rerun raw data or call fragments UMIs.
p='work/oct1-deep/neoantigen/bap1-local-validation/results.json';sources.add(p);support=[]
for r in j(p)['results']:
 for x in r['peptide_haplotype_support']:
  f=x['fragments'];rd=x['reads'];sd=x['read_strands'];support.append(dict(gene='BAP1',assay=r['kind'],scope=x['label'],alt_reads=rd.get('exact_deletion_haplotype',0),ref_reads=rd.get('exact_reference_haplotype',0),other_reads=rd.get('other_haplotype',0),alt_queryname_fragments=f.get('exact_deletion_haplotype',0),ref_queryname_fragments=f.get('exact_reference_haplotype',0),other_queryname_fragments=f.get('other_haplotype',0),discordant_fragments=f.get('discordant',0),alt_forward_reads=sd.get('exact_deletion_haplotype_forward',0),alt_reverse_reads=sd.get('exact_deletion_haplotype_reverse',0),source=p))
p='outputs/caris-deep-review/LATS-independent-read-audit.json';sources.add(p)
for r in j(p)['findings']:
 for kind,x in r['assays'].items():
  f=x['fragments'];rd=x['reads'];sd=x['strands'];support.append(dict(gene=r['gene'],assay=kind,scope='exact 17bp reference-span local haplotype',alt_reads=rd.get('alt',0),ref_reads=rd.get('ref',0),other_reads=rd.get('other',0),alt_queryname_fragments=f.get('alt',0),ref_queryname_fragments=f.get('ref',0),other_queryname_fragments=f.get('other',0),discordant_fragments=f.get('discordant',0),alt_forward_reads=sd.get('alt_forward',0),alt_reverse_reads=sd.get('alt_reverse',0),source=p))
tsv('read-support.tsv',support)
# All 31 equal-length near-self controls retained. Fixed allele, no new model inference.
ap='outputs/caris-deep-review/technical-evidence/12-current-reference-sensitivity/reference-control-affinity.tsv';np='outputs/caris-deep-review/technical-evidence/12-current-reference-sensitivity/gencode50-neighbors.tsv';sources.update([ap,np]);ar=list(csv.DictReader((ROOT/ap).open(),delimiter='\t'));nr={r['reference_peptide']:r for r in csv.DictReader((ROOT/np).open(),delimiter='\t')};controls=[]
for seq in dict.fromkeys(r['peptide'] for r in ar):
 scores={r['model_release']:r for r in ar if r['peptide']==seq};n=nr.get(seq);row=scores['2.3.0'];p20=scores['2.2.0'];controls.append(dict(peptide=seq,role=row['role'],allele='HLA-B*50:01',affinity_nM_model2_3_0=row['predicted_affinity_nM'],percentile_model2_3_0=row['affinity_percentile'],affinity_nM_model2_2_0=p20['predicted_affinity_nM'],percentile_model2_2_0=p20['affinity_percentile'],hamming_to_EERKGLYL=(sum(a!=b for a,b in zip(seq,'EERKGLYL')) if len(seq)==8 else ''),gencode50_occurrences=n['reference_occurrences'] if n else '',example_gencode50_header=n['example_header'] if n else '',sequence_reference_only='Yes' if n else 'No',experiment='paired dose-response specificity control' if n else 'candidate or same-position WT control'))
assert len(controls)==35 and len(nr)==31
fasta('BAP1-peptides-and-controls.fasta',[dict(id=f'peptide_{i+1:02d}',description=f'{r["role"]}; HLA-B*50:01 predicted only; peptide={r["peptide"]}',sequence=r['peptide']) for i,r in enumerate(controls)])
tsv('BAP1-peptides-and-controls.tsv',controls)
for src,dest in [('outputs/caris-deep-review/technical-evidence/03-HLA/patient-HLA-comparison.tsv','HLA-research-calls.tsv'),('outputs/caris-deep-review/Peptide-model-provenance.json','peptide-model-provenance.json'),('outputs/caris-deep-review/technical-evidence/02-BAP1-sequence-and-antigen/exon4-splice-junction-evidence.tsv','BAP1-splice-evidence.tsv')]:sources.add(src);shutil.copy2(ROOT/src,O/dest)
(O/'sequence-checks.json').write_text(json.dumps(dict(method='Independent local GenBank parsing, explicit strand-aware genomic edit and standard-code translation; conditional single-edit reference context',checks=checks),indent=2)+'\n')
manifest=[dict(path=p,sha256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest(),bytes=(ROOT/p).stat().st_size) for p in sorted(sources)]
(O/'source-manifest.json').write_text(json.dumps(dict(created='2026-10-04',note='No new patient uploads; same BAM analyses are computational replication, not laboratory validation',sources=manifest),indent=2)+'\n')
print(json.dumps(checks,indent=2));print('Wrote',len(list(O.iterdir())),'files')
