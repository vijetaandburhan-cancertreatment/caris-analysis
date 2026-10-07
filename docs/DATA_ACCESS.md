# Data access and privacy

This private repository contains patient-derived genomic information. Variant tables, read names, peptide sequences and selected BAMs are sensitive even without a name or a clinical report. Repository access is intended for authorized review of this case; it is not permission to publish, redistribute or upload the data to another analysis service.

## Original inputs

Nine Caris files total **33,619,145,590 bytes**. The three small originals are included in [data/small-inputs](../data/small-inputs/). Six large BAM/FASTQ files remain outside Git history, alongside the full set in the [controlled-access Drive folder](https://drive.google.com/drive/folders/1qQ_CQtDxoMX8TKYfPYFg0NSaKI2l900r). Individual file URLs, exact names, sizes and SHA256 values are in [raw-inputs.json](../manifests/raw-inputs.json) and [raw-inputs.tsv](../manifests/raw-inputs.tsv).

GitHub membership does not automatically grant Drive access. Request access from the repository owner if a linked input cannot be opened. Keep the original filenames and verify every downloaded file with `scripts/verify_inputs.py` before analysis. Do not replace an original with a newly processed file using the same name.

The repository is the analysis and provenance layer; controlled file storage is the large-data layer. Files beyond GitHub's size limits do not need to be split, truncated or recompressed merely to fit Git. The manifest connects both layers by byte-level checksums.

## What has actually been verified

- All nine resident original files were fully read and checked against their recorded SHA256 and source CRC values on 1 October 2026. All nine local SHA256 values were checked again successfully on 6 October 2026.
- All nine Drive copies had exact byte sizes checked.
- Only the three small Drive files received full download-and-SHA256 round-trip checks. The six large Drive copies were **not independently rehashed after upload** in that verification record.

This distinction matters: matching remote size alone is weaker evidence than a full checksum. A reviewer who downloads the large files can close that remaining transfer-verification gap by comparing their hashes with the manifest. Repository bundle checksums separately verify the packaged research files; they do not extend the scope of the historical Drive check.

## Sharing a review

Keep the GitHub repository and raw-data folder private unless the owner explicitly authorizes a different arrangement. Do not paste patient reads or derived results into public GitHub issues, public notebooks, public model endpoints or third-party web tools. Review findings can be returned in the private repository with the relevant commit, methods and evidence paths.

This repository contains no blanket grant to redistribute third-party software, references or model weights. Obtain those from the upstream sources and comply with their licenses. Historical records may retain local paths or accessions for provenance; removing a patient name alone would not anonymize the genomic data.
