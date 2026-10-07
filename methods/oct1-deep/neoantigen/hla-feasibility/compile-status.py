from pathlib import Path
import json,hashlib,time
OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[3]
BASE=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/hla-tools')
ARC=BASE/'arcasHLA';CTRL=BASE/'arcasHLA-control-3.24'
control=json.loads((OUT/'same-version-control.json').read_text());newer=json.loads((OUT/'public-genotype-smoke.json').read_text());smoke=json.loads((OUT/'extraction-smoke.json').read_text())
refs=[]
for label,base in [('3.46.0',ARC),('3.24.0_control',CTRL)]:
    for n in ['hla.idx','hla.fasta','hla.p.json']:
        p=base/'dat/ref'/n;refs.append({'version':label,'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
status={
    'created_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
    'status':'standard_genotyping_setup_and_extraction_mechanics_verified_no_patient_genotype',
    'arcas_source_commit':'9fa54a212d134b0d9894d1fc19ec1bdc6f62eb55',
    'arcas_scripts':str(ARC/'scripts'),
    'python':'/Users/burhanazeem/.local/share/codex/caris-analysis/bioenv/bin/python',
    'kallisto044':str(ROOT/'work/oct1-deep/fusion-next-step/public/kallisto-v0.44.0/kallisto/kallisto'),
    'requires_DYLD_override':False,
    'reference_files':refs,
    'total_task_tool_directory_bytes_including_small_smoke_outputs':sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file()),
    'resource_observations':{'reference_builder_python_peak_RSS_bytes':1239826432,'combined_python_kallisto_observed_RSS_approx_GB':1.7,'note':'Observed footprint, not an instrumented full process-tree peak. Optional partial index briefly raised generated reference directory to1.51GB, then only generated optional partial files were removed; final core/control resources below1GB. No original patient data removed or modified.'},
    'same_version_public_regression':control['status'],
    'same_version_seven_two_field_pairs_all_match':all(z['two_field_match'] for z in control['comparison'].values()),
    'newer_reference_public_control_caveat':'IMGT3.46 standard genotypes differ from the old3.24 complete-stage control at B,DQB1,DRB1; six of seven agree with upstream final partial-stage expectation, DRB1 differs. Do not hide reference sensitivity or call this a fully validated clinical pipeline.',
    'patient_mechanics_control':next(z for z in smoke['results'] if z['label']=='patient_HLA_A_mechanics_only'),
    'public_test_quality_caveat':'Upstream public BAM has no quality strings. Public-test-only placeholder qualities support sequence-only kallisto input; no fabricated patient qualities and no public base-quality inference.',
    'changes_to_official_algorithms':'Official downloaded source files unchanged. Reference-building wrappers replace only git-version lookup with explicit pinned direct-download commit and suppress optional partial-index output for resource limit. Extraction is a documented local pysam adapter, separately tested.',
    'next_bounded_step':'Whole intended RNA HLA subset (chr6 plus explicitly included unmapped-read sensitivity), proper pair/orphan accounting and input checks, then exploratory arcasHLA classI/classII with prior/no-prior review. Compare classI to Caris, retain classII ambiguity and reference-version sensitivity. This has not been run.',
    'limitations':['No full patient HLA extraction or genotype performed.','HLA-A region mechanics subset cannot support whole-panel genotype.','RNA genotype is not tumor allele-retention or LOH measurement; normal-cell admixture and allele silencing confound it.','Two-field inference remains provisional; neither a third field printed by software nor a strong gene TPM establishes HLA allele certainty.','Class-II inference may add a CD4/long-peptide research screen using a dedicated predictor; it does not establish tumor HLA-II expression, antigen presentation or efficacy.','FFPE/short reads, no matched normal, absent molecular purity/ploidy fit and older reference coverage limit clinical inference.'],
    'evidence_files':['feasibility.txt','input-feasibility.json','public-reference-manifest.json','optional-partial-cleanup.json','extraction-smoke.json','public-genotype-smoke.json','same-version-control.json'],
    'source_urls':['https://github.com/RabadanLab/arcasHLA','https://pmc.ncbi.nlm.nih.gov/articles/PMC6956775/','https://raw.githubusercontent.com/RabadanLab/arcasHLA/master/environment.yml','https://raw.githubusercontent.com/RabadanLab/arcasHLA/master/scripts/align.py','https://raw.githubusercontent.com/pachterlab/kallisto/v0.51.1/src/main.cpp','https://link.springer.com/article/10.1186/s13073-023-01154-x','https://www.nature.com/articles/nature22991','https://github.com/mskcc/lohhla','https://github.com/DiltheyLab/HLA-LA'],
    'script_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.glob('*.py')}
}
(OUT/'setup-status.json').write_text(json.dumps(status,indent=2)+'\n')
print(json.dumps({k:status[k] for k in ['status','same_version_public_regression','total_task_tool_directory_bytes_including_small_smoke_outputs']},indent=2))
