"""Conditional coding-sequence screen; no somatic status or immunogenicity claim."""
from pathlib import Path
import csv,gzip,hashlib,json,re,time,urllib.request,urllib.parse,io
from Bio import SeqIO
from Bio.Seq import Seq

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
OLD=ROOT/'work/oct1-analysis/transcripts'
SOURCES=list(csv.DictReader((ROOT/'outputs/caris-analysis/Caris-candidate-evidence.tsv').open(),delimiter='\t'))
REF=OUT/'reference-transcripts';REF.mkdir(exist_ok=True)
accessions=sorted({r['transcript'] for r in SOURCES})
missing=[]
for a in accessions:
    old=OLD/(a+'.gb');p=REF/(a+'.gb')
    if old.exists() and not p.exists():p.write_bytes(old.read_bytes())
    if not p.exists():missing.append(a)
if missing:
    url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode({'db':'nuccore','id':','.join(missing),'rettype':'gb','retmode':'text'})
    data=urllib.request.urlopen(url,timeout=90).read().decode()
    (REF/'retrieved-batch.gb').write_text(data)
    records=list(SeqIO.parse(io.StringIO(data),'genbank'))
    assert {r.id for r in records}==set(missing)
    for r in records:SeqIO.write(r,REF/(r.id+'.gb'),'genbank')

alleles=['HLA-A*02:01','HLA-A*26:01','HLA-B*50:01','HLA-B*52:01','HLA-C*12:02','HLA-C*06:02']
variants=[];peptides=[]
for src in SOURCES:
    record=SeqIO.read(REF/(src['transcript']+'.gb'),'genbank')
    assert record.id==src['transcript']
    cds_features=[f for f in record.features if f.type=='CDS'];assert len(cds_features)==1
    f=cds_features[0];cds=str(f.extract(record.seq));wt=str(Seq(cds).translate(cds=True))
    assert wt==f.qualifiers['translation'][0]
    snv=re.fullmatch(r'c\.(\d+)([ACGT])>([ACGT])',src['coding'])
    if snv:
        n,ref,alt=snv.groups();n=int(n);assert cds[n-1]==ref
        edited=cds[:n-1]+alt+cds[n:]
    else:
        spec=next(r for r in json.loads((OLD/'transcript-consequences.json').read_text())['results'] if r['gene']==src['gene'])
        a,b=spec['start_c'],spec['end_c'];assert cds[a-1:b]==spec['ref']
        edited=cds[:a-1]+spec['alt']+cds[b:]
    mutant=str(Seq(edited[:len(edited)//3*3]).translate(to_stop=True))
    first=next((i for i,(a,b) in enumerate(zip(wt,mutant)) if a!=b),min(len(wt),len(mutant)))
    pos=int(re.search(r'\d+',src['protein']).group());assert first+1==pos
    match=re.fullmatch(r'p\.([A-Z])(\d+)([A-Z*])',src['protein'])
    if match:
        ra,p,ma=match.groups();p=int(p)-1
        assert wt[p]==ra and (mutant[p] if p<len(mutant) else '*')==ma
    ident=f"{src['gene']}_{src['protein'].replace('*','Ter')}"
    v={'id':ident,'gene':src['gene'],'transcript':src['transcript'],'reported_protein':src['protein'],'coding':src['coding'],'wt':wt,'mutant':mutant,'first_changed_0based':first,'RNA_alt_fragments':int(src['strict_RNA_alt_fragments']),'somatic_status':'unknown; matched normal unavailable','input_record_sha256':hashlib.sha256((REF/(src['transcript']+'.gb')).read_bytes()).hexdigest()}
    variants.append(v)
    # Peptides must contain at least one non-reference amino acid. Stops have none.
    for length in range(8,12):
        for start in range(max(0,first-length+1),len(mutant)-length+1):
            end=start+length
            mt=mutant[start:end];refp=wt[start:end]
            if mt==refp:continue
            assert len(mt)==length and end>first
            peptides.append({'candidate_id':f'{ident}_{start+1}_{end}','variant':ident,'gene':src['gene'],'transcript':src['transcript'],'reported_protein':src['protein'],'start_1based':start+1,'end_1based':end,'length':length,'peptide':mt,'wildtype_corresponding_peptide':refp,'n_flank':mutant[max(0,start-5):start],'c_flank':mutant[end:end+5],'wt_n_flank':wt[max(0,start-5):start],'wt_c_flank':wt[end:end+5],'RNA_alt_fragments':v['RNA_alt_fragments'],'somatic_status':v['somatic_status']})

proteome=OUT/'gencode.v37.pc_translations.fa.gz'
queries={p['peptide']:[] for p in peptides}|{p['wildtype_corresponding_peptide']:[] for p in peptides}
proteins=0;residues=0
with gzip.open(proteome,'rt') as fh:
    for rec in SeqIO.parse(fh,'fasta'):
        proteins+=1;seq=str(rec.seq);residues+=len(seq)
        for q,hits in queries.items():
            # Exact string screen, all annotated isoforms. No approximate matching.
            i=seq.find(q)
            if i>=0:hits.append({'reference_id':rec.id,'first_match_start_1based':i+1})
for p in peptides:
    hits=queries[p['peptide']];p['normal_reference_exact_match_count']=len(hits)
    p['normal_reference_exact_match_ids']=';'.join(h['reference_id'] for h in hits)
    p['wildtype_reference_exact_match_count']=len(queries[p['wildtype_corresponding_peptide']])

with (OUT/'candidate-peptides.tsv').open('w') as out:
    w=csv.DictWriter(out,fieldnames=list(peptides[0]),delimiter='\t');w.writeheader();w.writerows(peptides)
fields=['candidate_id','variant','gene','sequence_type','allele','peptide','n_flank','c_flank','normal_reference_exact_match_count']
with (OUT/'mhcflurry-input.csv').open('w') as out:
    w=csv.DictWriter(out,fieldnames=fields);w.writeheader()
    for p in peptides:
        for typ in ['mutant','wildtype']:
            for allele in alleles+[';'.join(alleles)]:
                w.writerow({'candidate_id':p['candidate_id'],'variant':p['variant'],'gene':p['gene'],'sequence_type':typ,'allele':allele,'peptide':p['peptide'] if typ=='mutant' else p['wildtype_corresponding_peptide'],'n_flank':p['n_flank'] if typ=='mutant' else p['wt_n_flank'],'c_flank':p['c_flank'] if typ=='mutant' else p['wt_c_flank'],'normal_reference_exact_match_count':p['normal_reference_exact_match_count'] if typ=='mutant' else p['wildtype_reference_exact_match_count']})

manifest={'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'reference_url':'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.pc_translations.fa.gz','reference_sha256':hashlib.sha256(proteome.read_bytes()).hexdigest(),'reference_bytes':proteome.stat().st_size,'reference_translation_records':proteins,'reference_amino_acid_residues':residues,'input_variants':len(variants),'candidate_peptides':len(peptides),'unique_mutant_sequences':len({p['peptide'] for p in peptides}),'exact_normal_match_candidates':sum(p['normal_reference_exact_match_count']>0 for p in peptides),'zero_normal_match_candidates':sum(p['normal_reference_exact_match_count']==0 for p in peptides),'wildtype_controls_missing_from_reference':sum(p['wildtype_reference_exact_match_count']==0 for p in peptides),'per_variant':[{k:v[k] for k in ['id','gene','transcript','reported_protein','coding','RNA_alt_fragments','input_record_sha256']}|{'candidate_count':sum(p['variant']==v['id'] for p in peptides),'normal_match_count':sum(p['variant']==v['id'] and p['normal_reference_exact_match_count']>0 for p in peptides)} for v in variants],'limitations':['Sequences assume named RefSeq transcript, reference initiation and one stated variant; full-length transcript/phasing/protein not measured.','No matched-normal data: novel against reference is not evidence that an alteration is tumor-specific.','GENCODE37 does not cover every normal polymorphism, unannotated ORF, or processed peptide.','Exact absence does not establish safety, HLA binding, presentation, immunogenicity or clinical value.','Stop-gain variants have no new amino-acid sequence in this simple translation model and are not enumerated as novel peptides.','HLA alleles are tumor-inferred Caris results requiring blood confirmation.','Paired-read identifiers are not independent UMI molecules.']}
(OUT/'candidate-preparation.json').write_text(json.dumps(manifest,indent=2)+'\n')
(OUT/'conditional-proteins.json').write_text(json.dumps(variants,indent=2)+'\n')
(OUT/'normal-reference-matches.json').write_text(json.dumps({q:h for q,h in queries.items() if h},indent=2)+'\n')
print(json.dumps({k:v for k,v in manifest.items() if k!='per_variant'},indent=2),flush=True)
print(json.dumps(manifest['per_variant'],indent=2),flush=True)
