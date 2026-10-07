# Terms and measures

| Term | Meaning in this repository |
|---|---|
| Vendor call | A label or result supplied by Caris, such as “Deleted.” Repeating that label does not independently validate its biological interpretation. |
| Research observation | Something counted or reconstructed from the supplied reads, under stated filters and reference assumptions. |
| Model hypothesis | A prediction, such as potential HLA binding. A high score does not establish a working therapy. |
| MTAP / CDKN2A / CDKN2B | Three genes in the chromosome 9p21 region. The exported deletion labels do not by themselves reveal calibrated tumor copy numbers. |
| BAP1 / RASA1 / LATS1 / LATS2 | Genes with coding observations discussed here. A gene's known pathway role does not prove the functional effect of a particular observed change. |
| REF / ALT; WT / mutant | Reference/alternate nucleotide allele; reference (“wild type”)/changed sequence. Reference status is not the same as healthy-cell status. |
| VAF | Variant allele fraction: the proportion of qualifying sequence observations supporting a variant. Tumor purity, copy number and mixed cells affect it; it is not directly the fraction of cancer cells. |
| TPM | Transcripts per million, a normalized RNA abundance measure. Bulk TPM is not tumor-specific protein abundance. |
| HLA | Molecules involved in presenting peptides to immune cells. RNA-based HLA inference is provisional here; it does not establish tumor-cell HLA retention. |
| Peptide / 8mer / 9mer | A short amino-acid sequence; 8mer and 9mer indicate eight and nine amino acids. Sequences reconstructed from DNA/RNA are not necessarily produced or displayed by the tumor. |
| EL / BA | CapHLA's eluted-ligand-related and binding-affinity-related model scores. EL is not a patient-specific presentation probability; the normalized BA score is not a measured affinity in nM. |
| Binding register / core | The position of a candidate binding segment within a longer peptide. Possible cores listed here are not measured binding registers. |
| LOH / biallelic loss | Loss of heterozygosity / loss or inactivation of both gene copies. Tumor-only allele imbalance does not establish either mechanism by itself. |
| Query-name count / fragment | Paired reads grouped by their sequencing name. This avoids counting two mates as two fragments; it does not establish independent original molecules. |
| UMI | Unique molecular identifier: a molecular barcode. Query names in these counts are not UMIs. |
| MAPQ / BQ | Aligner mapping-quality and base-quality fields used in read filters. Their meanings depend on the software; an aligner's special value is not automatically a calibrated probability. |
| BAM / FASTQ / VCF | Alignment records / read sequences and qualities / variant records. These files contain patient-derived genomic data. |
| Fusion / splice junction | A joined sequence or RNA exon boundary. A supported junction need not be a genomic rearrangement, driver or treatment target. |
| SHA256 | A checksum used to compare exact file bytes. Agreement checks file integrity; it does not prove a scientific interpretation. |
| Internal cross-check | A separate implementation or review pass within the AI-assisted workflow. It is not outside peer review or independent laboratory validation. |

Read the [current findings](FINDINGS.md) with their adjoining limitations, and the [reproduction guide](REPRODUCTION.md) before running historical scripts.
