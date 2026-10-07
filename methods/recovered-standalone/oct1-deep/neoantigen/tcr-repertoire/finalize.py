"""Readable research findings; preserve executed binaries and output hashes."""
from pathlib import Path
import hashlib,json
OUT=Path(__file__).resolve().parent
BASE=Path('/Users/burhanazeem/.local/share/codex/caris-analysis/trust4-tools')
SRC=BASE/'TRUST4-a3fedd4aa0c1ad4da815427d82e0ceb69ce9c3a0'
x=json.loads((OUT/'findings.json').read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
x['executed_binary_hashes']=[{'name':name,'path':str(SRC/name),'bytes':(SRC/name).stat().st_size,'sha256':sha(SRC/name)} for name in ['trust4','bam-extractor','fastq-extractor','annotator']]
x['setup_note']='Initial clang build could not locate standard C++ vector header. A task-local explicit installed SDK/include configuration succeeded without source edits, toolchain installation or global system changes; both build logs/commands retained.'
x['unmapped_handling']='BAM is coordinate-sorted, but sampled unmapped records were not reliably in adjacent mate pairs. Used the official documented --abnormalUnmapFlag, which directs name-based candidate handling. No patient BAM modification.'
x['handoff_strength']='Limited. One TRB CDR3 comparison candidate has two source-verified Q30 fragment names. This is not a robust functional TCR discovery or a validated TIL result. Remaining reconstructed rows are exploratory and mostly single-fragment or otherwise insufficient under the explicit filter.'
(OUT/'findings.json').write_text(json.dumps(x,indent=2)+'\n')
lines=[
 'Exploratory bulk RNA T-cell receptor recovery — 1 October 2026','',
 'The existing bulk FFPE RNA contains a small amount of recoverable TCR sequence information. TRUST4 reported 17 TCR rows: 6 TRB, 7 TRA and4 TRG. These are reconstructed sequence candidates, not17 validated T-cell clones.','',
 'One TRB CDR3 passes the deliberately conservative local comparison filter: candidate TCR0002, annotated TRBV7-9*01 / TRBJ1-1*01. Its exact CDR3 nucleotide sequence is supported by2 different original paired-name fragments, both Q30 or better throughout the CDR3, with at least5 bases between the CDR3 and either read end. The2 full-read sequences differ. Both match the original BAM sequence and quality exactly. The CDR3 is unambiguous, in frame, stop-free and not imputed. The assembly is NOT complete VDJ. Evidence is limited to2 fragments; PCR-independent molecules have not been established.','',
 'The practical use is to retain this CDR3 nucleotide/amino-acid sequence plus V/J annotation for comparison with future paired single-cell TCR data from Eta/Michael. A future sequence match could identify a corresponding chain for further study. Neither the match nor this bulk sequence establishes antigen specificity, tumor reactivity or therapeutic usefulness. Do not manufacture or pair a receptor from this partial chain.','',
 'No TRA candidate reaches the two-source-verified-fragment threshold. The other reconstructed TCR rows remain in all-reconstructed-TCR.tsv, with their counts and limitations, so they can be rechecked if better data become available. They have not been silently discarded or promoted.','',
 'Quality and provenance: the official pinned TRUST4 source (v1.1.11-r641, commit a3fedd4aa0c1ad4da815427d82e0ceb69ce9c3a0) was built locally against installed Apple SDK headers, without source edits or global changes. The bundled public BAM control reproduces all44 CDR3 nucleotide identities, amino-acid labels, read counts and frequencies; exact full rows differ in2 consensus IDs and1 singleton IGHJ6 allele label. This is execution validation using a BCR-only example, not a TCR sensitivity benchmark. The official evaluation repository did not contain a small ready-to-run TCR input.','',
 'Patient run: unchanged RNA BAM; unchanged bundled hg38_bcrtcr.fa and human_IMGT+C.fa;1 thread; documented --abnormalUnmapFlag and --outputReadAssignment. All69,411 extracted read pairs have matching names and valid sequence/quality lengths, with no repeated pair names. Direct CDR3 matching collapsed mates by query name and counted original-quality Q20/Q30 evidence. All23 Q20-supporting mates across the reconstructed TCR candidates were found in indexed chr7/chr14/unmapped BAM scans and matched byte-for-byte in original orientation.','',
 'Resources: approximately4 minutes47 seconds for the patient TRUST4 run; sampled process-tree peak RSS311 MB. Additional local tools and outputs remained under100 MB. Original files preserved. No patient sequences were uploaded.','',
 'Interpretive limits:']+['- '+s for s in x['limitations']]+['',
 'Files: findings.json contains full reproducibility metadata and candidate details; read-supported-TRA-TRB.tsv contains the one limited comparison candidate; all-reconstructed-TCR.tsv preserves all17 TCR rows. Raw TRUST4 outputs, extracted reads, references and executed binaries remain under '+str(BASE)+'.','',
 'Primary method and implementation sources:',
 'https://github.com/liulab-dfci/TRUST4/tree/a3fedd4aa0c1ad4da815427d82e0ceb69ce9c3a0',
 'https://www.nature.com/articles/s41592-021-01142-2']
(OUT/'findings.txt').write_text('\n'.join(lines)+'\n')
print(json.dumps({'finalized':True,'files':['findings.json','findings.txt','all-reconstructed-TCR.tsv','read-supported-TRA-TRB.tsv'],'additional_directory_bytes':sum(p.stat().st_size for p in BASE.rglob('*') if p.is_file())},indent=2))
