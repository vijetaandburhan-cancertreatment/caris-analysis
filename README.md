# Caris DNA/RNA research — TN26-279853

Private evidence and reproducibility repository for independent review of nine supplied Caris bulk tumor DNA/RNA files. Analysis snapshot: **5 October 2026**; repository assembled **6 October 2026**.

**Start here:** [one-page findings (PDF)](reports/Caris-DNA-RNA-research-handoff-TN26-279853-2026-10-05-v2.pdf) · [plain text](reports/Caris-DNA-RNA-research-handoff-TN26-279853-2026-10-05-v2.txt) · [findings with evidence links](docs/FINDINGS.md)

The main starting points are MTAP/CDKN2A/CDKN2B copy-number interpretation; expressed BAP1/RASA1 coding events and provisional antigen hypotheses; additional LATS observations; RNA fusion/splice review; and DNA/RNA allele consistency under independent alignment. Observed read support, computational predictions and biological validation are explicitly distinguished. These analyses do not establish a diagnosis, prescribe treatment or predict clinical benefit.

## What is included

| Location | Contents |
|---|---|
| [`reports/`](reports/) | Latest five research summaries, plus clearly separated historical reports. |
| [`evidence/computational-followups-v2/`](evidence/computational-followups-v2/) | Exact allele/peptide packages, copy-number/phase work, HLA RNA and HLA-II investigations, fusion/splicing checks, controls and independent audits. |
| [`evidence/population-panel/`](evidence/population-panel/) | Frozen population-SNP panel, per-locus counts, annotations and original DNA/RNA comparison. |
| [`evidence/alignment-comparison/`](evidence/alignment-comparison/) | Fixed-cohort comparison with independently realigned RNA, retained read evidence and audits. |
| [`methods/`](methods/) | Additional methods and as-executed scripts, including recovered standalone discovery methods; more scripts are preserved inside evidence packages. |
| [`data/small-inputs/`](data/small-inputs/) | Original Caris VCF, gene-TPM CSV and analytical workbook, unchanged. |
| [`manifests/raw-inputs.json`](manifests/raw-inputs.json) | All nine original filenames, sizes, SHA-256 checksums and controlled-access Drive links. |
| [`docs/`](docs/) | Reproduction guidance, methods/versions, evidence limits and data-access instructions. |

## Large raw inputs

The nine originals total **33,619,145,590 bytes (33.62 GB)**. The six BAM/FASTQ files remain in the [private raw-data Drive folder](https://drive.google.com/drive/folders/1qQ_CQtDxoMX8TKYfPYFg0NSaKI2l900r); the three small originals are also included here. Local originals remain preserved separately.

This is a two-part research package: the versioned repository plus the exact raw inputs identified by its manifest. GitHub access does **not** automatically grant access to Drive. Reviewers need permission from the maintainers for both. No patient-derived files should be made public by changing repository or Drive visibility.

Each downloaded file should be verified against the manifest before analysis. All nine local originals were fully hashed and matched the source CRC on 1 October; all nine local SHA-256 checks passed again on 6 October. All nine Drive sizes were verified on 1 October. Only the three small Drive files received an independent full-download SHA-256 check; the large remote files were not independently rehashed.

## Quick verification

Python 3's standard library is sufficient for these integrity checks:

```sh
python3 scripts/verify_bundle.py
python3 scripts/verify_inputs.py --data-dir data/small-inputs --small-only
# After downloading the nine originals to a separate directory:
python3 scripts/verify_inputs.py --data-dir /path/to/raw-inputs
```

The first command verifies manifest-listed repository file bytes against recorded checksums. It does not check unlisted extra files, prove provenance or validate scientific claims. Do not edit historical evidence files when reproducing work; save new outputs separately and state the repository commit being reviewed.

## Reproducing and extending the work

Read [reproduction instructions](docs/REPRODUCTION.md) and [methods and versions](docs/METHODS_AND_VERSIONS.md). This is a collection of executed research workflows and their evidence, **not a tested portable one-command pipeline**. Original scripts preserve absolute paths and run assumptions. Configure paths and dependencies deliberately before rerunning.

The completed full-library genome-aware RNA run used STAR/Arriba with a sparse D8/SA12 index and the GENCODE v37 GRCh38 primary reference. It processed all 23,209,264 RNA pairs. A complete newly aligned patient BAM was not retained; selected read records, caller/junction outputs and audits were. The separate prepared D1 cloud workflow was not executed and is not represented as completed work.

Reports and evidence are versioned historical artifacts. Some include absolute local links or execution notes; the repository's relative links and manifests are the navigation layer. Earlier findings can be superseded by the 5 October v2 reports and later fixed-cohort alignment comparison. See [limitations](docs/LIMITATIONS.md) before interpreting a negative result or prioritizing a biological hypothesis.

The [method recovery record](manifests/method-recovery.json) supplements the earlier assembly record: it distinguishes methods already represented in the evidence packages, additional recovered scripts and intentional exclusions. The current repository file manifest covers the final combined package.

For an independent review, report the exact inputs, reference versions, commands, controls and differences from the current findings. Useful unresolved questions are listed at the end of [FINDINGS.md](docs/FINDINGS.md). Computational replication and biological or clinical validation are different outcomes.

See [data access](docs/DATA_ACCESS.md) and [reuse notice](NOTICE.md).
