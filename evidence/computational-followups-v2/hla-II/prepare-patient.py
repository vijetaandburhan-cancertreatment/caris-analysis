import pathlib,json,hashlib,itertools,pandas as pd
P=pathlib.Path(__file__).resolve().parent
src=P.parents[1]/'oct1-deep/neoantigen/BAP1-RASA1-local-protein-context.json'
data=json.loads(src.read_text()); lib=pd.read_csv(P/'CapHLA/HLA_library.csv').set_index('Allele Name')['MHC pseudo-seq'].to_dict()
alleles=[('HLA-DRB1*03:01','Provisional RNA beta-chain; DRA01:01 assumed'),('HLA-DRB1*15:02','Provisional RNA beta-chain; DRA01:01 assumed'),('HLA-DRB5*01:02','Sole RNA sequence type, not proven genomic dosage; DRA01:01 assumed')]
alleles += [(f'HLA-DPA1*01:03/DPB1*{b}','DPA1 sole RNA sequence type; DPB1 provisional RNA allele; not clinical typing') for b in ['02:01','04:01']]
alleles += [(f'HLA-DQA1*{a}/DQB1*{b}','Unphased exploratory alpha/beta combination; cis/trans or functional pairing not established') for a in ['01:03','05:01'] for b in ['02:01','06:01']]
assert all(a in lib for a,_ in alleles)
rows=[];windows=[]
for gene in ['BAP1','RASA1']:
 rec={r['sequence_type']:r for r in data['records'] if r['gene']==gene};m=rec['mutant'];w=rec['wt'];offset=m['first_altered_position1']-m['amino_acid_start1'];assert m['sequence'][:offset]==w['sequence'][:offset]
 for length in range(13,26):
  for start in range(len(m['sequence'])-length+1):
   end=start+length
   if end<=offset:continue
   assert m['sequence'][start:end]!=w['sequence'][start:end]
   wid=f'{gene}_{m["amino_acid_start1"]+start}_{length}'
   record={'window_id':wid,'gene':gene,'length':length,'aa_start1':m['amino_acid_start1']+start,'aa_end1':m['amino_acid_start1']+end-1,'altered_positions_in_peptide1':list(range(max(offset-start+1,1),length+1)), 'reaches_predicted_stop':end==len(m['sequence']),'mutant':m['sequence'][start:end],'wt':w['sequence'][start:end], 'transcript':m['transcript'],'coding_variant':m['coding_variant']}
   windows.append(record)
   for a,caveat in alleles:
    for kind in ['mutant','wt']:
     rows.append({'peptide':record[kind],'Allele Name':a,'window_id':wid,'gene':gene,'sequence_type':kind,'length':length,'aa_start1':record['aa_start1'],'aa_end1':record['aa_end1'],'altered_positions_in_peptide1':','.join(map(str,record['altered_positions_in_peptide1'])),'reaches_predicted_stop':record['reaches_predicted_stop'],'HLA_caveat':caveat})
pd.DataFrame(rows).to_csv(P/'patient-input.csv',index=False)
(P/'patient-peptide-manifest.json').write_text(json.dumps({'context_source':str(src),'context_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'context_records':data['records'],'alleles':[{'name':a,'pseudo_sequence':lib[a],'caveat':c} for a,c in alleles],'windows':windows,'total_windows':len(windows),'input_rows':len(rows),'scope':'Conditional mutation-containing13–25mers with same-start/same-length WT comparators. Do not extend natural new protein stops. DQ pairs unphased; no nearest HLA substitution.'},indent=2))
print(len(windows),len(rows),pd.DataFrame(windows).groupby('gene').size().to_dict())
