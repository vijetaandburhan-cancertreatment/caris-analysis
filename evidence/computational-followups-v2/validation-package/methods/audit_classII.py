from pathlib import Path
import csv,json,hashlib
P=Path('work/oct4-followup/hla-II');rows=list(csv.DictReader((P/'patient-paired-comparisons.csv').open()));manifest=json.loads((P/'patient-peptide-manifest.json').read_text());windows={x['window_id']:x for x in manifest['windows']}
assert len(rows)==1053 and len(windows)==117
assert len({(r['window_id'],r['Allele Name']) for r in rows})==1053
assert sum(w['gene']=='BAP1' for w in windows.values())==91
for r in rows:
 w=windows[r['window_id']];assert r['peptide']==w['mutant'] and r['WT_peptide']==w['wt']
 altered=[i for i,(a,b) in enumerate(zip(w['mutant'],w['wt'])) if a!=b]
 cores=[w['mutant'][s:s+9] for s in range(len(w['mutant'])-8) if any(s<=i<s+9 for i in altered)]
 assert len(cores)==int(r['mutation_containing_9mer_windows']);assert len(w['mutant'])-8-len(cores)==int(r['wholly_WT_9mer_windows'])
 assert cores==r['mutation_containing_9mer_sequences'].split(';')
lead=next(r for r in rows if r['window_id']=='RASA1_237_15' and r['Allele Name']=='HLA-DPA1*01:03/DPB1*02:01')
assert lead['peptide']=='GDYYIGGRRFSSLQT' and lead['WT_peptide']=='GDYYIGGRRFSSLSD'
assert max(float(r['v2_mutant_el']) for r in rows if r['gene']=='BAP1')<.371
result=dict(scope='Independent candidate row/count and possible-core arithmetic audit; not a model rerun or biological validation',paired_comparisons=1053,unique_mutant_windows=117,all_row_sequences_and_core_counts_match=True,lead=lead,source_sha256=hashlib.sha256((P/'patient-paired-comparisons.csv').read_bytes()).hexdigest())
(P/'validation-package-independent-audit.json').write_text(json.dumps(result,indent=2)+'\n');print('All 1053 comparison rows, 117 sequences and possible-core counts verified; RASA1 lead confirmed.')
