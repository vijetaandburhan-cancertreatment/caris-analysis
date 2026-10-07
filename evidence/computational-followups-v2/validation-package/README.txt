BAP1 / LATS research validation package | 4 October 2026

Start with BAP1-LATS-validation-handoff.pdf (four pages) or the identical plain-text narrative.
Machine-readable TSV/FASTA files carry exact versioned sequence definitions, controls and scoped read evidence. Full interpretation limitations are in the narrative.

These are reference-derived, single-edit research sequences, not clinically validated primers, fully phased patient transcripts, synthesis constructs, proven antigens or treatment instructions. All three variants need matched-normal/orthogonal confirmation.

methods/ contains the package-generation code, separate Biopython sequence audit and PDF renderer. Generation refers to pinned pre-existing workspace sources identified in source-manifest.json. The small relevant audit/reference sources are included in evidence/ and references/. Large raw sequence files are intentionally outside this handoff. No new patient upload or outreach occurred.

The base handoff contains no validated class-II peptide claim. supplements/ adds the controlled exploratory RASA1 CD4 lead with its paired reference control, and independently audited local BAP1 SNP/deletion phase. These addenda do not establish peptide presentation, TCR specificity, LOH or biallelic inactivation.

Integrity: checksums.sha256 covers every package file except itself. Check from inside this directory with: shasum -a 256 -c checksums.sha256
