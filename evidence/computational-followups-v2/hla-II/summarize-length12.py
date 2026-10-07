import pandas as pd,json,pathlib
P=pathlib.Path(__file__).resolve().parent;r=json.loads((P/'normal-reference-exact-check.json').read_text());m=json.loads((P/'length12-peptide-manifest.json').read_text());nov={s for s,x in r['cores'].items() if x['mutation_containing'] and x['reference_occurrences']==0};assert all(any(s in w['mutant'] for s in nov) for w in m['windows'])
frames=[]
for version in ['v2','v1']:
 d=pd.read_csv(P/f'length12-{version}-scores.csv');d['model_release']=version;frames.append(d)
f=pd.concat(frames);f.to_csv(P/'length12-all-scores.csv',index=False);out={'scope':'Post-primary sensitivity adding12mers because they are within the new method\'s published classII length range. Original13–25mer primary results preserved. No selection of favorable rows.','windows':9,'BAP1_windows':7,'RASA1_windows':2,'HLA_combinations':9,'rows_per_release':162,'reference_absence_inference':'Each new12mer contains a mutation-containing9mer with zero exact matches in the completedGENCODE50 scan; consequently no complete exact12mer reference match can exist. This is a logical consequence of the earlier scan, not a new full scan.','models':{}}
for version,d in f.groupby('model_release'):
 a=d[d.sequence_type=='mutant'].set_index(['window_id','Allele Name']);b=d[d.sequence_type=='wt'].set_index(['window_id','Allele Name']);a=a.copy();a['WT_peptide']=b.peptide;a['WT_EL']=b.presentation_score;a['WT_BA']=b.affinity_score
 out['models'][version]={'max_EL_by_gene':a.groupby('gene').presentation_score.max().to_dict(),'max_BA_by_gene':a.groupby('gene').affinity_score.max().to_dict(),'both_author_cutpoints_by_gene':a.assign(both=(a.presentation_score>.5)&(a.affinity_score>.5)).groupby('gene').both.sum().astype(int).to_dict(),'top_EL_per_gene':a.reset_index().sort_values('presentation_score',ascending=False).groupby('gene').head(1).to_dict('records')}
(P/'length12-sensitivity.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
