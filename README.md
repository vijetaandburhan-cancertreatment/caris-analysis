# Caris DNA/RNA research — TN26-279853

Public research methods, findings and genomic data from nine supplied Caris bulk tumor DNA/RNA files. Analysis snapshot: **5 October 2026**; repository assembled and annotated **6 October 2026**. The work is shared to enable scrutiny, replication and new analyses. It does not establish a diagnosis, prescribe treatment or predict clinical benefit.

**Start with [current findings and evidence](docs/FINDINGS.md).** For a brief overview, see the [one-page handoff](reports/Caris-DNA-RNA-research-handoff-TN26-279853-2026-10-05-v2.pdf) alongside the [current corrections and status map](docs/CURRENT_STATUS_AND_ERRATA.md). [Terms and measures](docs/GLOSSARY.md) explains the main terminology.

## What the work found

The main starting points are the vendor MTAP/CDKN2A/CDKN2B deletion calls and their unresolved copy-number interpretation; expressed BAP1/RASA1 coding events and provisional antigen hypotheses; additional unconfirmed LATS observations; RNA fusion/splice review; and DNA/RNA allele consistency under a second alignment pipeline. Observed reads, model predictions and biological validation are different evidence levels. The [findings table](docs/FINDINGS.md) states the limitations beside each result.

## How the work was produced

This is an **AI-assisted, family-led research analysis**. AI agents helped write and execute code, interpret outputs and prepare the reports. Checks described in historical files as “independent” used separate implementations within the same AI-assisted workflow. They are internal computational cross-checks, not external peer review, independent biological samples or laboratory validation. External expert review is invited and has not been established by this repository. See the [disclosure and interpretation limits](docs/LIMITATIONS.md).

## Review the evidence

| Location | Contents |
|---|---|
| [reports/](reports/) | Five current dated summaries and a separate historical handoff; read the current errata alongside them. |
| [evidence/computational-followups-v2/](evidence/computational-followups-v2/) | Allele/peptide packages, copy-number/phase work, HLA RNA and HLA-II work, fusion/splicing checks and controls. |
| [evidence/population-panel/](evidence/population-panel/) | Frozen public-SNP panel, per-locus patient counts, annotations and original DNA/RNA comparison. |
| [evidence/alignment-comparison/](evidence/alignment-comparison/) | Fixed-cohort comparison using the same RNA reads under another alignment pipeline, with retained records and internal audits. |
| [evidence/full-caller-outputs/](evidence/full-caller-outputs/) | Actual STAR/Arriba run record, caller/junction outputs and selected BAM records. |
| [methods/](methods/) | As-executed scripts and provenance; further scripts are inside the evidence packages. |
| [data/small-inputs/](data/small-inputs/) | Unchanged original Caris VCF, gene-TPM CSV and analytical workbook. |
| [manifests/raw-inputs.json](manifests/raw-inputs.json) | All nine original filenames, sizes, SHA256 values and download links. |

## Get the full raw data

All nine original files total **33,619,145,590 bytes (33.62 GB)**. The complete set is linked from the [public raw-data Drive folder](https://drive.google.com/drive/folders/1qQ_CQtDxoMX8TKYfPYFg0NSaKI2l900r). Three small originals are also in this repository; six large BAM/FASTQ files remain outside Git because of file-size limits. Local originals remain preserved separately. This release includes patient-derived reads, genotypes and HLA-related results; it is not anonymized by the use of an accession identifier.

Download files using the [manifest](manifests/raw-inputs.json), retain their names and verify their bytes before analysis. [Data access and verification](docs/DATA_ACCESS.md) explains which checks were performed. Public accessibility, correct remote size and a range/header check do not establish the SHA256 of an entire large remote file.

```sh
python3 -m unittest discover -s tests -v
python3 scripts/check_release.py
python3 scripts/verify_bundle.py
python3 scripts/verify_inputs.py --data-dir data/small-inputs --small-only
# After downloading the nine originals:
python3 scripts/verify_inputs.py --data-dir /path/to/raw-inputs
```

These checks use the Python 3.9+ standard library. The 19 synthetic/negative tests check arithmetic, navigation and parser behavior; the release check compares 11 aggregate metrics and local Markdown link destinations. The bundle verifier rejects unlisted files as well as missing or mismatched manifest-listed files. The raw-input verifier checks the selected original files, without assessing unrelated directory contents. These checks do not rerun the genomic pipeline, establish scientific correctness or replace biological validation.

## Reproduce or extend the analysis

Read [reproduction instructions](docs/REPRODUCTION.md) and [recorded methods and versions](docs/METHODS_AND_VERSIONS.md). This is a collection of executed workflows and their evidence, **not a tested portable one-command pipeline**. Original scripts retain absolute paths and run assumptions. A clean-machine end-to-end rerun has not been demonstrated by publishing this repository.

The full-library RNA run processed all **23,209,264 RNA pairs** using STAR/Arriba, a sparse D8/SA12 index and the GENCODE v37 GRCh38 primary reference. A complete newly aligned patient BAM was not retained; selected records and caller/junction outputs were. The prepared D1 cloud workflow was not executed. Numerical cross-checks do not remove these limits.

Dated reports and evidence remain unchanged historical records; some contain local links, earlier statuses or terminology clarified by the [current status and errata](docs/CURRENT_STATUS_AND_ERRATA.md). The [method recovery record](manifests/method-recovery.json) identifies additional recovered source and intentional exclusions. The current repository manifest describes the assembled files.

For a review, report the repository commit, input/reference hashes, commands, controls, numerical differences and unresolved interpretations. See [useful next validation work](docs/FINDINGS.md#useful-next-work-for-an-external-validator) and the [reuse notice](NOTICE.md).
