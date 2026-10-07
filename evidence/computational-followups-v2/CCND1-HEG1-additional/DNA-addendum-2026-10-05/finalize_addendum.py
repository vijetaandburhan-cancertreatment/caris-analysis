#!/usr/bin/env python3
import pathlib,json,hashlib,shutil,datetime,sys
P=pathlib.Path(__file__).resolve().parent;R=P.parent;M=pathlib.Path.home()/'.local/share/codex/caris-analysis/oct5-CCND1-HEG1-audit-v1'/P.name
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
now=lambda:datetime.datetime.now(datetime.timezone.utc).isoformat()
prior_expected='ef400711acf72ded0c49f4d2282d9359fe88d8abeb56fdce4447c6adf4d2471b'
assert sha(R/'checksums.json')==prior_expected
prior=json.loads((R/'checksums.json').read_text());prior_resident=M.parent;prior_mismatches=[n for n,h in prior.items()if sha(prior_resident/n)!=h];assert not prior_mismatches,prior_mismatches
for n in('independent-root-DNA.json','independent-linked-root-rescore-cases.json','independent-root-patient-only-rescore-cases.json','independent-DNA-linked-sequence-membership.json'):assert(P/n).is_file()
prov={'completed_utc':now(),'task':'Separate bounded indexed DNA addendum to completed CCND1::HEG1 RNA audit','prior_RNA_manifest_sha256':prior_expected,'prior_RNA_manifest_entries_verified':len(prior),'prior_RNA_reverification_source':str(prior_resident),'workspace_prior_RNA_manifest_unchanged':True,'original_RNA_results_unchanged':True,'original_DNA_full_file_not_rescanned_or_copied':True,'original_DNA_verification':'Size and mtime unchanged; source index SHA recorded in dna-audit.json. A new full-source-BAM SHA was not calculated for this bounded check.','only_large_reference_scan':'Two linked DNA sequence families plus public/reference-derived controls searched against existing local primary GRCh38 and G37 transcripts; no full patient DNA scan.','reuse_script_hashes':{n:{'original':sha(R/n),'addendum':sha(P/'whole-reference-alternatives'/n),'identical':sha(R/n)==sha(P/'whole-reference-alternatives'/n)}for n in('nominate.py','score.py')},'reference_sources':'Same preserved genomic FASTA and G37 transcript reference as prior RNA audit; local only.','independent_review':'Root independent marker/CIGAR/counts, 288 synthetic marker-BQ cases, 12 local score cases and 20 whole-reference top score cases; BAM membership checked separately.','clinical_validation':False,'no_external_contact_or_new_upload_in_this_addendum':True,'limits':'Exact marker negative plus two nonduplicate-name non-exact DNA hints. No confirmed genomic fusion, matched-normal status or driver/treatment claim.'}
(P/'provenance.json').write_text(json.dumps(prov,indent=2)+'\n')
exclude={'checksums.json','combined-receipt.json','mirror-verification.json'}
files=sorted(x for x in P.rglob('*')if x.is_file()and'__pycache__'not in x.parts and x.name not in exclude and x.suffix!='.pyc')
man={str(f.relative_to(P)):sha(f)for f in files};(P/'checksums.json').write_text(json.dumps(man,indent=2)+'\n')
receipt={'completed_utc':now(),'prior_RNA_version_manifest_sha256':prior_expected,'prior_RNA_entries_reverified':len(prior),'prior_RNA_entries_reverified_from_resident_mirror':True,'DNA_addendum_manifest_sha256':sha(P/'checksums.json'),'DNA_addendum_entries':len(man),'DNA_addendum_findings_sha256':sha(P/'findings.txt'),'DNA_addendum_bounded_BAM_sha256':sha(P/'bounded-DNA-records.bam'),'local_directory':str(P),'resident_mirror':str(M),'separate_versions':True,'manifest_self_hash_excluded':True}
(P/'combined-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
M.mkdir(parents=True,exist_ok=True)
copyfiles=files+[P/'checksums.json',P/'combined-receipt.json']
for f in copyfiles:
 dst=M/f.relative_to(P);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(f,dst)
assert all(sha(f)==sha(M/f.relative_to(P))for f in copyfiles)
verify={'verified_utc':now(),'files_compared':len(copyfiles),'bytes':sum(f.stat().st_size for f in copyfiles),'all_hashes_match':True,'mirror':str(M),'prior_RNA_manifest_unchanged':sha(R/'checksums.json')==prior_expected,'prior_RNA_manifest_in_resident_mirror_matches':sha(M.parent/'checksums.json')==prior_expected,'DNA_manifest_sha256':sha(P/'checksums.json'),'combined_receipt_sha256':sha(P/'combined-receipt.json')}
(P/'mirror-verification.json').write_text(json.dumps(verify,indent=2)+'\n');shutil.copyfile(P/'mirror-verification.json',M/'mirror-verification.json');assert sha(P/'mirror-verification.json')==sha(M/'mirror-verification.json')
print(json.dumps({'receipt':receipt,'verification':verify},indent=2))
