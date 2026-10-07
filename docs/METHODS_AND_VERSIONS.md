# Recorded methods, versions and references

These are versions recorded in the completed analyses, not a new unified environment specification. Exact commands and file hashes in the linked JSON records are authoritative for the corresponding run. A separate [upstream version check](../manifests/upstream-version-check.json) records official package/release sources and the model-asset hash check; this verifies version provenance, not biological accuracy. No new computational re-analysis or clean-machine end-to-end rerun is implied by packaging them here. Historical “independent” labels refer to internal cross-checks; see [disclosure and errata](CURRENT_STATUS_AND_ERRATA.md).

## Whole-reference RNA alignment and fusion review

| Component | Recorded execution | Source / provenance |
|---|---|---|
| STAR | **2.7.11b**, upstream macOS x86_64 binary, run on Apple silicon through Rosetta | [Native tool record](../methods/genome-fusion/native-tool-provenance.json), [actual patient command](../evidence/full-caller-outputs/run.json); [upstream](https://github.com/alexdobin/STAR) |
| Arriba | **2.5.1**, built natively for arm64 from upstream source | Same native record; [upstream](https://github.com/suhrig/arriba) |
| Genome | GRCh38 primary assembly distributed with **GENCODE v37** | [Reference provenance](../evidence/computational-followups-v2/fusion-robustness/reference-provenance-addendum.json) |
| Annotation | `gencode.v37.primary_assembly.annotation.gtf` | GENCODE v37 primary-assembly annotation, not a latest-version substitution |
| STAR index | **D8/SA12**, sparse index; six threads in the recorded patient fusion run | [Run command](../evidence/full-caller-outputs/run.json) and index-building source under [methods/genome-fusion](../methods/genome-fusion/) |
| Arriba resources | v2.5.1 GRCh38 blacklist, known-fusion recovery/tags and protein domains; default caller filters retained | [Execution provenance](../evidence/computational-followups-v2/fusion-robustness/execution-provenance-addendum.json) |

The Arriba build needed macOS SDK/include configuration and a uniquely named libdeflate ARM CPU-feature object to avoid a build-object basename collision. The provenance record states that algorithm source was not changed. A separate diagnostic build added counters; its patches, checks and comparisons are preserved under [caller-diagnostics](../evidence/computational-followups-v2/caller-diagnostics/).

The completed fusion run processed 23,209,264 RNA pairs: 81.06% uniquely mapped and 13.21% multilocus-mapped in STAR's reported metrics. It produced 28 accepted fusion rows and 825,436 discarded rows; 26 discarded hypotheses were selected under documented rules for deeper review. The run streamed all alignments; no full realigned BAM was saved.

The original vendor RNA BAM records STAR **2.7.8a**, a GRCh38 / GENCODE37 CTAT reference and two-pass settings. The local pipeline used STAR 2.7.11b, primary-genome D8 and one-pass settings, with additional overlap, indel, splice and chimeric differences. Original contig names/lengths agree with the local reference, but the exact original reference/index sequence identity is unavailable. See [original-versus-new provenance](../evidence/alignment-comparison/independent-comparison/original-new-STAR-provenance.json). This comparison cannot isolate sparsity or prove that the newer mapping is correct.

## Public reference downloads

| File | Official source | SHA256 of compressed download |
|---|---|---|
| GRCh38 primary genome | [GENCODE release 37 FASTA](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.primary_assembly.genome.fa.gz) | `142c225c317be446bd2005ce717db02cc7276fa6bc667d6fe423bd1e0fe64c99` |
| GENCODE v37 primary annotation | [GENCODE release 37 GTF](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.primary_assembly.annotation.gtf.gz) | `fa9e291d9e6e0b68db8c225239878f77c315f4a959f4fe81b42535d67eb9af56` |

The decompressed genome SHA256 recorded in the reference addendum is `e49b92b3e4f321bf254c042f25b726d9931c4d74c7523e8b6bb530e63b0cfd4b`. The download records also checked official MD5 values. Do not compare a compressed-file hash with its decompressed form.

The public SNP panel uses **1000 Genomes 2019 GRCh38**, biallelic PASS variants with global ALT frequency 0.05–0.95, on chromosomes 2, 3, 5, 8, 9 and 17. Population selection, intervals, complete accounting and exact frozen design are under [population-panel](../evidence/population-panel/); every selected REF was checked against the resident FASTA. Selection did not use this patient's DNA allele calls. SAS frequency was recorded, not used for selection.

## Base counting and quality controls

The population analysis records **samtools 1.24 / pysam 0.24.1**. Main filters: MAPQ ≥30, BQ ≥25, proper pairs, exclude flags `0xF0C`, BAQ disabled, samtools mate-overlap adjustment disabled, and explicit `-d 0` depth setting (qualified as effectively unlimited in the installed version). Qualifying A/C/G/T observations collapse by exact query name; conflicting mate bases are discarded. Reference skips/deletions are not A/C/G/T coverage. A clean-read sensitivity adds no soft clipping and five-base aligned-end margins; MAPQ60-only checks remain separate.

Installed depth behavior was checked with 9,001 synthetic paired names / 18,002 records and finite-cap contrasts. A separate internal CIGAR-counting implementation reproduced the 4,015 shared positions. The later alignment comparison used matched counting criteria with a fixed original cohort and retained coverage failures. See [population controls and counts](../evidence/population-panel/) and [comparison controls and audits](../evidence/alignment-comparison/).

## Peptide and HLA model work

| Component | Recorded version / provenance | Interpretation limit |
|---|---|---|
| MHCflurry | **2.3.9**; model release **2.3.0**; PyTorch 2.14.1 recorded. [Model provenance and hashes](../evidence/computational-followups-v2/validation-package/peptide-model-provenance.json); [upstream](https://github.com/openvax/mhcflurry) | Binding estimates are not demonstrated tumor presentation or recognition. |
| CapHLA v1 | Commit `a17016fae96eae81bd710384fa288c17e84d2d3d` | Public numerical controls reproduced; per-allele training counts unavailable. |
| CapHLA v2 | Commit `33ebdd6ce6dadbbb1c66b026ce4b5d81dbf3a831` | EL output is not a patient-specific probability. |
| CapHLA execution | Python **3.12.14**, PyTorch **2.14.1**, NumPy **2.3.5**, pandas **2.2.3** recorded; local CPU wrapper | [Method validation](../evidence/computational-followups-v2/hla-II/method-validation.json), [version record](../evidence/computational-followups-v2/hla-II/model-version-provenance.json); [upstream](https://github.com/changyunjian/CapHLA) |

CapHLA v1's supplied allele table lacked pseudosequences expected by its code. The official v2 mapping was used without inventing residues, and all 100 released v1 score rows, including 26 DP rows, reproduced within the recorded tolerance. Architectures/utilities matched across versions; weights were separately identified. HLA class-II inference remains exploratory: some chains are assumed, DQ pairings unphased, and one model pseudosequence retains an upstream `X`.

An earlier DeepSeqPanII regression did **not** meet its checks; no patient scores were produced with it. The failed check is part of method-selection history, not validated antigen evidence. The replacement was selected because public regression passed, not because of a favorable patient result.

Third-party model weights, full reference databases, tool environments and prebuilt indexes are not generally vendored here. Their provenance allows a reviewer to obtain the specified versions. Historical methods and supplemental source are under [methods](../methods/); check each workstream's manifest rather than assuming all scripts ran in one environment.

## Technical benchmark scope

The public BCR::ABL1 experiment comprised **45 technical runs**, including 20 matched compact D1/D8 conditions representing 19 unique paired-input sets. It was not a 45-patient benchmark. The small correlated mixtures, known-fusion priors and limited background cannot estimate clinical sensitivity or a detection limit. See the [benchmark evidence](../evidence/computational-followups-v2/fusion-robustness/) and [current errata](CURRENT_STATUS_AND_ERRATA.md).
