"""Read only the first 10,000 pairs; preserve all original FASTQs and prior QC."""
from pathlib import Path
import collections
import datetime
import gzip
import json
import re
import numpy as np

source = Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
out = Path(__file__).parent
names = [f'RNA_TN26-279853_S25.R{i}.fastq.gz' for i in (1, 2)]
start = datetime.datetime.now(datetime.timezone.utc).isoformat()
lengths = [collections.Counter(), collections.Counter()]
metadata = [collections.Counter(), collections.Counter()]
qualities = [[], []]
sequences = [[], []]
first_headers = [[], []]
barcode_tags = [collections.Counter(), collections.Counter()]
mate_labels = [collections.Counter(), collections.Counter()]
format_errors = []
mate_mismatches = 0
pairs = 0
eof_reached = False
files = [gzip.open(source / n, 'rb') for n in names]
try:
    for k in range(10000):
        records = [[f.readline().rstrip(b'\r\n') for _ in range(4)] for f in files]
        if not records[0][0] or not records[1][0]:
            eof_reached = True
            if bool(records[0][0]) != bool(records[1][0]):
                format_errors.append({'pair': k + 1, 'error': 'Asymmetric EOF'})
            break
        ids = []
        for j, (header, sequence, plus, quality) in enumerate(records):
            if not header.startswith(b'@') or not plus.startswith(b'+') or len(sequence) != len(quality):
                format_errors.append({'pair': k + 1, 'mate': j + 1, 'error': 'FASTQ framing or sequence/quality length mismatch'})
            parts = header.split(maxsplit=1)
            ids.append(parts[0].removesuffix(b'/1').removesuffix(b'/2'))
            suffix = parts[1].decode('ascii', 'replace') if len(parts) > 1 else ''
            metadata[j][suffix] += 1
            mate_labels[j][suffix.split(':', 1)[0] if suffix else 'absent'] += 1
            for tag in re.findall(rb'(?:^|\s)((?:CB|CR|UB|UR|RX|BC):[A-Za-z]:\S+)', header):
                barcode_tags[j][tag.decode('ascii', 'replace')] += 1
            if k < 2:
                first_headers[j].append(header.decode('ascii', 'replace'))
            lengths[j][len(sequence)] += 1
            qualities[j].append(quality)
            sequences[j].append(sequence)
        mate_mismatches += ids[0] != ids[1]
        pairs += 1
finally:
    for f in files:
        f.close()

mates = []
for j, name in enumerate(names):
    q_ascii = np.frombuffer(b''.join(qualities[j]), dtype=np.uint8)
    q = q_ascii.astype(np.int16) - 33
    bases = b''.join(sequences[j]).upper()
    mates.append({
        'file': name,
        'records_sampled': sum(lengths[j].values()),
        'length_counts': dict(sorted(lengths[j].items())),
        'length_summary': {
            'minimum': min(lengths[j]),
            'maximum': max(lengths[j]),
            'mean': sum(n*c for n,c in lengths[j].items()) / sum(lengths[j].values()),
            'records_shorter_than_20_nt': sum(c for n,c in lengths[j].items() if n < 20),
            'records_shorter_than_30_nt': sum(c for n,c in lengths[j].items() if n < 30),
        },
        'bases_sampled': len(bases),
        'base_counts': {b: bases.count(b.encode()) for b in 'ACGTN'},
        'other_base_count': len(bases) - sum(bases.count(b.encode()) for b in 'ACGTN'),
        'quality_encoding_assumption': 'Phred+33, consistent with Illumina headers; not an independent assay metadata confirmation',
        'quality_ascii_min': int(q_ascii.min()),
        'quality_ascii_max': int(q_ascii.max()),
        'mean_phred_Q': float(q.mean()),
        'fraction_bases_Q20_or_higher': float((q >= 20).mean()),
        'fraction_bases_Q30_or_higher': float((q >= 30).mean()),
        'phred_outside_0_to_93_count': int(((q < 0) | (q > 93)).sum()),
        'first_headers': first_headers[j],
        'header_metadata_counts': dict(metadata[j]),
        'leading_mate_metadata_token_counts': dict(mate_labels[j]),
        'explicit_CB_CR_UB_UR_RX_BC_header_tags': dict(barcode_tags[j]),
    })

result = {
    'started_at': start,
    'finished_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'source_directory': str(source),
    'sampling': 'First 10,000 pairs sequentially; not a random or whole-file sample',
    'pairs_sampled': pairs,
    'EOF_reached_in_sample': eof_reached,
    'pair_ID_mismatches': mate_mismatches,
    'format_errors': format_errors,
    'sample_format_and_pairing_pass': pairs == 10000 and not mate_mismatches and not format_errors,
    'mates': mates,
    'barcode_interpretation': [
        'Read headers use matching instrument/flowcell/tile coordinates and distinct leading mate tokens 1 and 2.',
        'Metadata suffixes are preserved exactly; the nonstandard 1.5b component is not interpreted as a molecular or cellular barcode.',
        'No cell barcode or UMI identity can be established from matched read IDs or these metadata fields alone.',
        'No explicit CB/CR/UB/UR/RX/BC tags were found if the reported tag dictionaries are empty. Inline barcode/UMI sequence cannot be ruled out without the library protocol.',
        'These paired FASTQs and the existing STAR whole-transcriptome BAM header are consistent with bulk RNA sequencing. They do not establish single-cell barcodes or paired single-cell TCR data.'
    ],
    'processing_provenance_followup': 'Variable read lengths, a few one-base reads and the nonstandard metadata suffix justify clarifying pre-export trimming/processing. This is not proof of corruption or sample mismatch, and hashes match Caris source objects.',
    'full_raw_FASTQ_QC_status': 'Pending',
    'why_full_QC_pending': [
        'Source CRC64NVME plus local SHA256 cover every compressed byte and establish faithful transfer; they do not parse every decompressed FASTQ record, pairing, read length or quality.',
        'This limited preflight does not reach gzip EOF and therefore does not independently perform gzip trailer validation of the complete decompressed stream.',
        'Full BAM scans verify the delivered alignment representation, which may reorder, filter, clip or otherwise transform raw records; they are not a raw FASTQ all-record audit.',
        'Unmeasured full-file metrics include total FASTQ record count, pairing/order across the full files, full length/quality/per-cycle/adapter distributions and raw-to-BAM count reconciliation.',
        'The requested scope is a 10,000-pair preflight now. A later compiled-tool streaming full-file audit can compute those metrics without creating decompressed copies.'
    ],
    'originals_modified': False,
    'implementation': 'Python record loop; NumPy vectorized quality statistics and native bytes counting, no Python per-base loop'
}
(out / 'fastq-rna-preflight.json').write_text(json.dumps(result, indent=2) + '\n')
lines = [
    'RNA paired FASTQ preflight (first 10,000 pairs only)',
    f'Source: {source}',
    f'Pairs sampled: {pairs}; mate ID mismatches: {mate_mismatches}; format errors: {len(format_errors)}',
]
for m in mates:
    s = m['length_summary']
    lines.append(f"{m['file']}: length range {s['minimum']}-{s['maximum']} nt, mean {s['mean']:.2f} nt, {s['records_shorter_than_20_nt']} reads <20 nt; mean Q {m['mean_phred_Q']:.2f}; Q30 {100*m['fraction_bases_Q30_or_higher']:.2f}%; N bases {m['base_counts']['N']}; explicit barcode/UMI header tags {m['explicit_CB_CR_UB_UR_RX_BC_header_tags']}")
lines += ['', 'Barcode/read architecture:', *result['barcode_interpretation'], result['processing_provenance_followup'], '', 'Why full raw FASTQ QC remains pending:', *result['why_full_QC_pending']]
(out / 'fastq-rna-preflight.txt').write_text('\n'.join(lines) + '\n')
print('\n'.join(lines[:5]))
