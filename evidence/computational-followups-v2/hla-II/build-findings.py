import pathlib,json,pandas as pd,hashlib,datetime,shutil
P=pathlib.Path(__file__).resolve().parent; d=pd.read_csv(P/'patient-paired-comparisons.csv');s=json.loads((P/'summary.json').read_text());mv=json.loads((P/'method-validation.json').read_text())
r=d[(d.window_id=='RASA1_237_15')&(d['Allele Name']=='HLA-DPA1*01:03/DPB1*02:01')].iloc[0].to_dict()
q={}
for v in ['v1','v2']:
 x=pd.read_csv(P/f'patient-{v}-scores.csv');q[v]={'BA_min':float(x.affinity_score.min()),'BA_max':float(x.affinity_score.max()),'BA_outside_target_0_1':int(((x.affinity_score<0)|(x.affinity_score>1)).sum()),'execution':json.loads((P/f'patient-{v}-scores.csv.json').read_text())}
res={'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'Completed controlled exploratory HLA-II screen; no biological validation','method_validation':{'v2_public_examples':32,'v1_public_examples':100,'v1_classII_DP_examples':26,'max_absolute_tolerance':1e-5,'passed':True,'detail':'method-validation.json'},'scope':{'mutant_long_peptides':117,'BAP1':91,'RASA1':26,'lengths_aa':[13,25],'HLA_combinations':9,'mutant_WT_pairings_per_release':1053,'scored_rows_per_release':2106,'releases':2,'natural_stops_preserved':True},'gene_summary':s['gene_summary'],'research_lead':r,'regression_range_diagnostics':q,'normal_reference':s['normal_reference'],'prioritization':'RASA1 15mer/DP02:01 is an exploratory CD4 testing lead; high WT scores and shared WT cores prevent a mutation-specific claim. No robust BAP1 classII result. Does not invalidate earlier separate classI BAP1 predictions.','limits':['RNA HLA calls are provisional; DRA assumed; DQ combinations unphased; sole DPA1/DRB5 calls do not establish homozygosity.','EL scores are not calibrated patient presentation probabilities; BA is unbounded normalized regression, not measured nM or percentile rank.','The two releases are the same model family, not independent validation. Public numerical reproduction establishes execution fidelity, not clinical accuracy.','Allele-specific training counts could not be verified; v1 bundled classII controls are DP only.','Single-edit reference transcript translation; no matched normal or protein, HLA ligand or T-cell assay confirmation.','Mutated9mer windows are possible cores only; no binding register inferred. Exact reference absence does not establish tumor specificity or safety.','No patient data sent to external prediction service; all inference blocked network access.']}
res['length12_sensitivity']=json.loads((P/'length12-sensitivity.json').read_text())
(P/'findings.json').write_text(json.dumps(res,indent=2))
text='''HLA-II/CD4 follow-up — completed 4 October 2026

The most useful new result is a conditional RASA1 CD4 testing lead, not a validated neoantigen. BAP1 does not gain a robust class-II result from this screen; its earlier, separate class-I findings are unchanged.

We established a reproducible local replacement for the failed DeepSeqPanII run: CapHLA, using its official MIT-licensed source and models. All 32 released v2 test examples and all 100 v1 examples, including 26 class-II DP examples, reproduced within the prespecified 0.00001 tolerance (maximum error about 0.00000041). The earlier DeepSeqPanII mismatch remains documented. V1's released HLA table omitted sequences; supplying the official v2 sequence mapping reproduced its published predictions without modifying model weights or architecture.

The completed scan covers every mutation-containing 13–25mer possible in the existing BAP1/RASA1 contexts: 91 BAP1 and 26 RASA1 peptides, each with a matching normal-reference peptide, across nine provisional HLA-II combinations. That is 2,106 input rows per release and 1,053 mutant/WT comparisons. Natural stops are retained: seven altered amino acids for BAP1, two for RASA1. No downstream residues were invented.

RASA1: the natural-terminal 15mer GDYYIGGRRFSSLQT (amino acids 237–251), paired with provisional DPA1*01:03/DPB1*02:01, has a consistently high model EL score:
  Current v2: mutant 0.925750; WT GDYYIGGRRFSSLSD 0.777399.
  Older v1: mutant 0.958507; WT 0.802943.
  Normalized BA: v2 0.453212 versus WT 0.400247; v1 0.428404 versus WT 0.387869.
This is suitable for discussing a paired mutant-versus-WT CD4 assay with the research team after HLA confirmation. It is not evidence that T cells recognize the mutant selectively. Five of its seven possible nine-amino-acid binding cores are entirely WT; only GGRRFSSLQ and GRRFSSLQT contain the changes, at the ends of the possible cores. The actual bound register is unknown. Neither model establishes natural tumor presentation.

BAP1: within the primary 13–25mer scan, the maximum current-release EL score is 0.370006. Two older-release peptide/HLA pairs cross the authors' exploratory EL/BA 0.5 cutpoints, but both fail those combined cutpoints in the current release. No BAP1 or RASA1 pair passes both cutpoints in both releases. These thresholds are descriptive model screening rules, not a biological exclusion test; their failure does not prove that a peptide cannot work.

A separately labeled 12mer sensitivity added seven BAP1 and two RASA1 windows (162 rows per release). It found no pair passing both author cutpoints. BAP1 IFLFKWIEERKG / DRB5*01:02 scores EL 0.532382 in v2, but its WT IFLFKWIEERRS is almost identical at 0.511394, and BA is 0.273520 versus WT 0.279429. The older release has a larger EL contrast, making the apparent mutant selectivity unstable. This does not change the priority above. Across primary and sensitivity work, 126 distinct mutant windows of lengths 12–25 were examined, with matched WT controls.

Normal-reference check: none of the 117 mutant long peptides or nine mutation-containing possible nine-amino-acid cores exactly matches any of 382,428 GENCODE 50 translation records. This does not account for the patient's constitutional variants, unannotated proteins, alternate processing or structural cross-reactivity, and does not establish safety.

Keep the score meanings narrow. EL is a model output, not a percentage likelihood of presentation or benefit. BA is an unbounded normalized regression: 4/2,106 current-release and 45/2,106 older-release outputs fall outside the nominal 0–1 target interval. We retained them without clamping or converting them into apparently precise affinities. Specific training counts for Papa's alleles were not available. RNA HLA calls need clinical confirmation; DRA is assumed, all four DQ combinations are unphased, and two DQ mappings contain the upstream model's X residue. There is no matched-normal, protein, HLA-ligand or functional T-cell validation. All patient scoring was local, with network access blocked during inference.

Files: patient-paired-comparisons.csv contains all comparisons and fold ranges; per-HLA-candidate-comparison.csv shows the strongest EL, strongest BA and largest EL mutant/WT contrast for every gene/HLA pair. method-validation.json, the scripts, source manifests and README.txt make the run reproducible.

Primary method: https://github.com/changyunjian/CapHLA
Publication: https://doi.org/10.1093/bib/bbae595
'''
(P/'findings.txt').write_text(text)
selected=d[(d.window_id=='RASA1_237_15')&d['Allele Name'].isin(['HLA-DPA1*01:03/DPB1*02:01','HLA-DPA1*01:03/DPB1*04:01'])].copy();selected.to_csv(P/'RASA1-CD4-research-lead.csv',index=False)
# Compact deliverable; all large downloaded public sources and original outputs retained in work.
out=P.parents[1].parent/'outputs/caris-followup/hla-II';out.mkdir(parents=True,exist_ok=True)
for f in ['findings.txt','findings.json','README.txt','method-validation.json','patient-paired-comparisons.csv','per-HLA-candidate-comparison.csv','RASA1-CD4-research-lead.csv','patient-peptide-manifest.json','normal-reference-exact-check.json','CapHLA-source-manifest.json','CapHLA-v1-source-manifest.json','model-version-provenance.json','run_caphla.py','prepare-patient.py','summarize.py','check-reference.py','download-models.py','build-findings.py','independent-sequence-audit.txt','length12-sensitivity.json','length12-all-scores.csv','length12-peptide-manifest.json','prepare-length12-sensitivity.py','summarize-length12.py','patient-input.csv','patient-v2-scores.csv','patient-v1-scores.csv','length12-input.csv','public-v1-expected.csv','public-v1-replay.csv','public-v2-replay.csv']:
 shutil.copy2(P/f,out/f)
files=[{'path':str(f.relative_to(out)),'size':f.stat().st_size,'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sorted(out.iterdir()) if f.is_file() and f.name!='manifest.json']
(out/'manifest.json').write_text(json.dumps({'created_utc':res['finished_utc'],'source_work_directory':str(P),'files':files},indent=2));print(out)
