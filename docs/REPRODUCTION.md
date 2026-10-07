# How to audit and reproduce this work

The repository preserves evidence, reports, as-executed scripts and provenance. It is **not a tested one-command pipeline**. Historical scripts contain absolute paths, directory dependencies, local environments and assumptions about prepared references. Do not run them unchanged expecting a portable workflow.

## 1. Verify the copy you received

Clone the private repository, record its commit ID, then run these standard-library Python checks from the repository root:

```sh
python3 scripts/verify_bundle.py
python3 scripts/verify_inputs.py --data-dir data/small-inputs --small-only
# Once all nine original inputs have been downloaded:
python3 scripts/verify_inputs.py --data-dir /absolute/path/to/the/nine/Caris/files
```

The first checks packaged repository files against [the repository manifest](../manifests/repository-files.sha256.json). The input checks compare the bundled small files or the downloaded full set against [raw-inputs.json](../manifests/raw-inputs.json). These checks establish byte identity to the recorded versions, not biological correctness. Full-input hashing reads approximately 34 GB and can take time. Obtain the six large inputs separately as described in [DATA_ACCESS.md](DATA_ACCESS.md); the three small originals are also in [data/small-inputs](../data/small-inputs/).

## 2. Audit the results without reprocessing every read

Begin with [FINDINGS.md](FINDINGS.md). The three evidence directories preserve the most recent completed analyses:

| Directory | What to inspect first | Questions the package can answer |
|---|---|---|
| [computational-followups-v2](../evidence/computational-followups-v2/) | Topic findings, exact sequence/variant tables, scripts, independent audit receipts and model controls | Are the reported counts, translated windows, local phase and candidate-selection decisions supported by the captured evidence? |
| [population-panel](../evidence/population-panel/) | `design-frozen.json`, `all-public-loci-counts.tsv`, `primary-summary.json`, `joint-callable.tsv`, `coverage-attrition.tsv`, exception and audit directories | Are site selection, denominator, read filters, overlap handling and discrepancies reproduced? |
| [alignment-comparison](../evidence/alignment-comparison/) | `design-frozen.json`, `baseline-fixed-4015.tsv`, `comparison/fixed4015-comparison.tsv`, `comparison/summary.json`, `independent-comparison/` | What changes when the same reads are aligned by the other recorded pipeline, retaining the original cohort and coverage losses? |

Use numerical tables and audit receipts before relying on PDF prose. Original reports may contain local absolute links or operational notes from the original work session; those links are historical provenance, not portable download instructions. Current repository paths and the raw-input manifest are the navigation/access layer.

## 3. Re-run a bounded analysis from original inputs

1. Choose a specific claim and freeze the expected inputs, coordinates, transcript/reference version, filtering and output metrics before looking at alternative results.
2. Verify the raw files and obtain the exact public references/model versions in [METHODS_AND_VERSIONS.md](METHODS_AND_VERSIONS.md) and the linked manifests. Reference/model downloads are separate from this repository.
3. Inspect the relevant scripts and their imports/path constants. Create a working copy with paths adapted to your environment, or independently implement the method. Record the patch and new commands; preserve the originals as the execution record.
4. Recreate the recorded environment for that workstream. There is no single verified dependency lock covering all historical scripts. Native macOS binaries, Rosetta execution, CPU wrappers and model-library repairs require explicit attention if moving to another platform.
5. Run the provided public/synthetic numerical controls before patient inference. Reproduce exact integer counts where the method is deterministic; use the prespecified numerical tolerance for model scores rather than claiming byte-for-byte model reproducibility across hardware.
6. Write new outputs to a separate location. Compare both successful and failed/depth-lost sites with the frozen baseline, retaining changes in denominator and selection. Record tool versions, reference hashes, commands, exit codes, checksums and unresolved differences.

The original scripts are evidence of what was run, not a promise that every dependency or intermediate reference is vendored. If a required external dependency cannot be reconstructed, document the gap rather than substituting a new method silently.

## 4. Repeat full-reference RNA work

The completed fusion analysis processed **23,209,264 paired RNA reads** with STAR 2.7.11b, the GRCh38 primary-assembly / GENCODE v37 reference and a sparse **D8/SA12** index, streaming the full alignment into Arriba 2.5.1. The [recorded run](../evidence/full-caller-outputs/run.json) preserves the actual command arrays. All recorded process exit codes were zero. The source entry points are under [methods/genome-fusion](../methods/genome-fusion/).

A complete realigned BAM was **not saved**. The original run retained a small four-gene BAM for targeted splicing, along with splice/chimeric/caller outputs. The later alignment comparison performed another full-library pass and selectively captured the fixed SNP cohort and all emitted records for preselected discrepancy names. Its ~64 MB BAM is useful for that audit, not a substitute for a full BAM.

The earlier population addendum correctly reported that the four-gene BAM was inadequate to evaluate its broad panel. The later [full-library comparison](../reports/Caris-final-alignment-comparison-2026-10-05.txt) addresses a larger fixed cohort. Do not count the earlier non-evaluable check as a successful replication.

A prepared full-reference dense D1 cloud pipeline was not executed and is not presented as completed patient work. Compact D1/D8 public-control comparisons did run; these are distinct from an unperformed dense full-patient analysis. A new dense-index or different-aligner run would be an extension, requiring its own reference, parameters, resources and controls.

## 5. Return a useful validation result

Include the repository commit, input/reference hashes, software/environment record, code or patch, exact commands, controls and a compact expected-versus-observed table. State whether your result is an exact replication, a parameter sensitivity test, an alternative interpretation or new independent biological evidence. For new findings, preserve supporting and contradicting reads and explain why an ordinary-reference/repeat/WT explanation is insufficient.

See [LIMITATIONS.md](LIMITATIONS.md) for inference boundaries and [DATA_ACCESS.md](DATA_ACCESS.md) before transferring patient-derived content to another service.
