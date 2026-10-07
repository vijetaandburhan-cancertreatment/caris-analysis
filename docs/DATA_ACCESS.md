# Raw data, access and verification

This is a full public genomic research release. It includes patient-derived variants, HLA-related results, read evidence and the original DNA/RNA inputs. These materials remain potentially identifying; removing names would not anonymize the genetic data.

## Download the nine original files

The originals total **33,619,145,590 bytes**. All nine are linked from the [public raw-data Drive folder](https://drive.google.com/drive/folders/1qQ_CQtDxoMX8TKYfPYFg0NSaKI2l900r). The three small originals are also included in [data/small-inputs](../data/small-inputs/). Six large BAM/FASTQ files remain outside Git history. Individual links, filenames, byte sizes and SHA256 values are in [raw-inputs.json](../manifests/raw-inputs.json) and [raw-inputs.tsv](../manifests/raw-inputs.tsv).

Google may display a warning that a large file cannot be virus-scanned; use its **Download anyway** action to retrieve the original file. This expected size-related warning is separate from the checksum verification below.

Keep the original filenames. Check downloaded files from the repository root:

```sh
python3 scripts/verify_inputs.py --data-dir /absolute/path/to/raw-inputs
```

A full verification reads approximately 34 GB. A browser success page, matching file size or successful range download does not replace this checksum check. If access fails or a checksum differs, report the affected filename and result to the maintainers; do not silently use the file as verified input.

The repository provides versioned methods/results; Drive holds large data that exceeds GitHub's file limits. The checksum manifest connects those two parts. No file needs to be truncated or split to preserve the original data.

## Verification scope

| Check | What the record establishes |
|---|---|
| Original local verification, 1 October 2026 | All nine resident originals were fully read, SHA256 checked and matched against recorded source CRC64NVME values. |
| Local recheck, 6 October 2026 | All nine original SHA256 values matched again. |
| Drive upload verification, 1 October 2026 | Exact byte sizes matched for all nine. The three small files received complete download-and-SHA256 round-trip checks. |
| Public-access check, 6 October New York / 7 October UTC | All nine files were accessible without authentication. The three small files passed complete anonymous-download SHA256 checks. For the six large files, HTTP 206 ranges reported the expected total sizes and the first/last 65,536 bytes matched local originals. See [public-raw-access.json](../manifests/public-raw-access.json). These range checks do not verify every byte of a large file. |
| Six large remote files | A complete remote re-download-and-SHA256 verification was not established by the earlier upload/access checks. Each reviewer should verify their downloaded copies against the manifest. |
| Repository bundle checks | `scripts/verify_bundle.py` checks the inventory and bytes against the manifest, rejecting missing, mismatched and unlisted files. It does not prove scientific correctness. |

The newly aligned **full RNA BAM was never retained**. Selected alignment BAMs and the original vendor BAMs are provided; they must not be mistaken for a complete local realignment. Public references, model weights and tool environments are not all vendored. Their recorded versions and upstream sources are in [METHODS_AND_VERSIONS.md](METHODS_AND_VERSIONS.md).

## Submit a useful review

Include the repository commit, exact input hashes, software/reference versions, commands, controls and expected-versus-observed results. A code replication, parameter sensitivity test and independent biological validation are different contributions. Unresolved read origins and failed controls should remain visible.

See [NOTICE.md](../NOTICE.md) for reuse and AI-assistance disclosure. Public release does not create a blanket license for third-party software, references or model weights. Historical files may retain original paths/accessions and earlier access instructions; use this page and the current manifests for the present release.
