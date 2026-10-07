"""Independent public data/model provenance inspection. No model inference or patient peptides."""
import csv,json,pathlib,math,hashlib,collections,datetime
root=pathlib.Path(__file__).resolve().parent
pub=root.parent/'class-II-feasibility/public'
manifest=json.loads((pub.parent/'public-source-manifest.json').read_text())
source_hashes=[]
for r in manifest['files']:
 p=pub.parent/r['path'];h=hashlib.sha256(p.read_bytes()).hexdigest()
 source_hashes.append({'path':r['path'],'matches_manifest':h==r['sha256'],'sha256':h,'bytes':p.stat().st_size})
assert all(r['matches_manifest']for r in source_hashes)
train=list(csv.DictReader((pub/'code_and_dataset/dataset/BD_2013_DATAPROVIDER_READY.txt').open(),delimiter='\t'))
def auc(rows):
 cutoff=1-math.log(500)/math.log(50000)
 pos=[float(r['pred'])for r in rows if float(r['real'])>cutoff]
 neg=[float(r['pred'])for r in rows if float(r['real'])<=cutoff]
 # Direct pairwise calculation avoids depending on ML metric libraries.
 value=sum((p>n)+0.5*(p==n)for p in pos for n in neg)/(len(pos)*len(neg))if pos and neg else None
 return {'IC50_binding_cutoff_nM':500,'normalized_cutoff':cutoff,'positives':len(pos),'negatives':len(neg),'ROC_AUC_from_published_predictions':value}
results=[]
for key,label in [('DRA01:01-DRB103:01','DRA*01:01-DRB1*03:01'),('DRA01:01-DRB115:02','DRA*01:01-DRB1*15:02')]:
 file=pub/'Models/BD2016_LOMO'/key/'test_result.txt'
 rows=list(csv.DictReader(file.open(),delimiter='\t'));tr=[r for r in train if r['HLA']==label]
 results.append({'HLA':label,'published_test_rows':len(rows),'unique_test_peptides':len({r['pep']for r in rows}),'test_peptide_lengths':dict(sorted(collections.Counter(len(r['pep'])for r in rows).items())),'BD2013_training_rows':len(tr),'BD2013_unique_training_peptides':len({r['sequence']for r in tr}),'BD2013_peptide_lengths':dict(sorted(collections.Counter(len(r['sequence'])for r in tr).items())),'test_peptides_also_in_BD2013_same_allele':len({r['pep']for r in rows}&{r['sequence']for r in tr}),**auc(rows)})
out={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'repo_commit':manifest['commit'],'source_file_hash_verification':source_hashes,'source_file_count':len(source_hashes),'source_file_bytes':sum(r['bytes']for r in source_hashes),'total_BD2013_training_rows':len(train),'HLA_specific_data':results,'interpretation':'These published LOMO test predictions assess allele-held-out behavior. They are not an independent evaluation of the different BD2013 general checkpoint used for patient hypotheses, and test overlap with BD2013 is explicitly counted. Source-file prediction reproduction by the running agent tests computational consistency, not biological validity.'}
(root/'public-data-audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k!='source_file_hash_verification'},indent=2))
