HLA-II/CD4 research follow-up — local, reproducible, unvalidated biological predictions

Scope
This branch replaces an unresolved DeepSeqPanII numerical-reproduction failure with an independently published open-source local tool, CapHLA. It does not repair or erase that earlier failure. No patient sequence was sent to a prediction service; a Python audit hook prevents socket operations during inference. These are research prioritization scores, not medical advice, measured binding, validated neoantigens, or a vaccine design.

Source and execution
CapHLA current v2 commit33ebdd6ce6dadbbb1c66b026ce4b5d81dbf3a831; older v1 commita17016fae96eae81bd710384fa288c17e84d2d3d. Official repository MIT license retained. Source manifests record URLs, SHA256 hashes and sizes. Each release uses its original5 EL and5 BA model weights and the unchanged upstream model architectures and inference functions. CPU strict checkpoint load, eval mode and original one-hot encoding/X padding. Our wrapper changes only data plumbing: no multiprocessing, bounded batch32, two CPU threads. It avoids the upstream unused pickle/rank loading and keeps all five individual scores.

The original v1 HLA CSV omitted pseudosequences despite its script requiring them. The official v2 mapping was therefore supplied explicitly, with no reconstructed residues. All100 v1 released public predictions, including26 class-II DP rows, match within the predeclared1e-5 maximum absolute tolerance. All32 v2 public examples (classI/mouse, not classII) also pass. These numerical checks establish faithful local execution on public examples; they do not establish accuracy on Papa's DR/DQ alleles. Model weights differ between releases; v1 is sensitivity, not independent biological validation. Specific training counts by patient allele are unavailable in the checked repositories. The current repository supplies all9 exact requested identifiers; no nearest-allele substitution was used.

Peptides and HLA
prepare-patient.py enumerates all mutation-containing13–25mers possible within the existing conditional BAP1 and RASA1 protein contexts, with same-start/same-length WT comparators. The natural new stop is retained: BAP1 has7 altered residues, RASA1 has2. No residues, linkers, or tails are invented. Result:91 BAP1 plus26 RASA1 windows, each scored as mutant andWT across9 provisional HLA-II combinations (2106 rows/release). Exact HLA assumptions and pseudosequences are in patient-peptide-manifest.json. RNA HLA typing is not accredited germline typing. DRA01:01 is assumed; DRB5 andDPA1 sole RNA sequence types are not proven homozygosity. All4 DQ alpha/beta combinations are explicitly unphased; functional/cis/trans pairing is unproved. The upstream DQA1*05:01 pseudosequence has1X, retained unchanged and flagged.

Score meanings
EL is the released model's softmax presentation score, not a calibrated probability that this patient's tumor displays the peptide. BA is the raw normalized regression score, higher predicting stronger binding; it is not nM or percentile rank. Scores are not vaccine efficacy, T-cell recognition, clinical benefit, or safety probabilities. Five-fold minima/maxima are descriptive model variation, not confidence intervals. Author cutpoints EL>0.5 andBA>0.5 are descriptive exploratory screen cutpoints only. Never relabel a high-scoring mutation-containing long peptide as a mutation-specific antigen:85/91 BAP1 windows and all26 RASA1 windows contain at least one wholly WT9mer. The algorithm does not give a validated binding register.

Normal-reference screen
check-reference.py searches exact longer peptides and their possible9mer cores against pinned GENCODE50 protein-coding translations (all382428 records, SHA2565d6408c3a1c22d864c96a55367cf5b69f1031dea1fe2cd6d03e366ae3fc56a45). No exact match for117 mutant long peptides or9 mutation-containing9mers. This is a reference translation collection, not a measured normal proteome or matched-normal genome. No absence-of-self or safety inference follows.

Reproduction (run from workspace root; use Python with versions in score JSON)
python work/oct4-followup/hla-II/download-models.py
python work/oct4-followup/hla-II/run_caphla.py --input work/oct4-followup/hla-II/CapHLA/test_out.csv --output work/oct4-followup/hla-II/public-v2-replay.csv --expected work/oct4-followup/hla-II/CapHLA/test_out.csv
python work/oct4-followup/hla-II/run_caphla.py --input work/oct4-followup/hla-II/public-v1-expected.csv --output work/oct4-followup/hla-II/public-v1-replay.csv --expected work/oct4-followup/hla-II/public-v1-expected.csv --model-dir work/oct4-followup/hla-II/CapHLA-v1 --library work/oct4-followup/hla-II/CapHLA/HLA_library.csv
python work/oct4-followup/hla-II/prepare-patient.py
python work/oct4-followup/hla-II/run_caphla.py --input work/oct4-followup/hla-II/patient-input.csv --output work/oct4-followup/hla-II/patient-v2-scores.csv
python work/oct4-followup/hla-II/run_caphla.py --input work/oct4-followup/hla-II/patient-input.csv --output work/oct4-followup/hla-II/patient-v1-scores.csv --model-dir work/oct4-followup/hla-II/CapHLA-v1 --library work/oct4-followup/hla-II/CapHLA/HLA_library.csv
python work/oct4-followup/hla-II/check-reference.py
python work/oct4-followup/hla-II/summarize.py

Primary sources
https://github.com/changyunjian/CapHLA
https://doi.org/10.1093/bib/bbae595
https://www.gencodegenes.org/human/
Additional limits
The code accepts7–25 amino acids; published length-specific class-II evaluation highlighted12–20mers. The21–25mer results are retained, but should not be promoted merely for a high regression score. Some affinity regressions can fall outside0–1; outputs are retained without clamping or conversion into apparent measured affinities. The peptide sequence contexts assume the selected reference transcript and single reported edit; they do not establish complete personalized cis haplotypes, full-length protein production or its abundance. Potential class-II priming through antigen-presenting cells is distinct from malignant-cell class-II expression, neither of which is demonstrated by these scores.

Compact delivery folder use
The outputs/caris-followup/hla-II folder contains the exact input tables and scores, plus code and pinned public manifests. To reproduce predictions without regenerating peptide contexts, change into that folder, install the package versions recorded in findings.json, run python download-models.py, then:
python run_caphla.py --input patient-input.csv --output reproduced-v2.csv
python run_caphla.py --input patient-input.csv --output reproduced-v1.csv --model-dir CapHLA-v1 --library CapHLA/HLA_library.csv
The download script retrieves public model files only. Input sequences remain on the machine. The full original generation and reference-check workflow uses the preserved work directory and source contexts; source records are also embedded in patient-peptide-manifest.json. The separate length12-input.csv can be scored with the same commands. That post-primary extension adds nine windows and is reported separately in length12-sensitivity.json, preserving the original13–25mer analysis.
