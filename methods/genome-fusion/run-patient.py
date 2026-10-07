"""All-pair RNA genome-aware exploratory analysis after positive public controls."""
import pathlib,sys,json,datetime
P=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(P));from align import run
C=pathlib.Path('/Users/burhanazeem/.local/share/codex/caris-analysis');R=C/'genome-fusion-reference';D=C/'TN26-279853'
controls=json.loads((P/'full-reference-controls/results.json').read_text());assert controls['status']=='complete' and controls['BCR_ABL1_execution_control_pass']
compact=json.loads((P/'compact-controls/results.json').read_text());assert compact['comparison']['D1_positive'] and compact['comparison']['D8_positive']
gate={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'public_execution_controls_passed':True,'sparse_equivalence_gate':'NOT PASSED: preserved compact D1/D8 differences. Proceed for discovery with explicit sensitivity limitation; independent full D1 comparison planned by parent.','filter_tuning':False,'all_expected_pairs':23209264,'notes':['Genome-aware RNA fusion analysis, not whole-genome DNA sequencing.','Arriba does not report intragenic deletions; negative results cannot exclude BAP1/LATS exon loss.','No full patient BAM retained locally; four-gene selective BAM and full junction/SJ logs retained.']}
(P/'patient-start-gate.json').write_text(json.dumps(gate,indent=2))
result=run(R/'STAR_G37_primary_D8_SA12',R/'GRCh38.primary_assembly.genome.fa',R/'gencode.v37.primary_assembly.annotation.gtf',D/'RNA_TN26-279853_S25.R1.fastq.gz',D/'RNA_TN26-279853_S25.R2.fastq.gz',P/'patient-D8',threads=6,save_bam=False,expected_pairs=23209264,capture_bed=P/'targeted-splicing/genes-GENCODE37.bed')
print(json.dumps({'status':result['status'],'finished_utc':result['finished_utc'],'STAR_final':result['STAR_final']},indent=2))
