"""Bounded, read-only RNA BAM fusion-input audit; no reference downloads."""
from pathlib import Path
import collections
import datetime
import json
import shutil
import subprocess
import pysam

out = Path(__file__).parent
source = Path('/Users/burhanazeem/.local/share/codex/caris-analysis/TN26-279853')
counts = collections.Counter()
tags = collections.Counter()
with pysam.AlignmentFile(str(source/'RNA_TN26-279853.bam'), 'rb') as bam:
    header = bam.header.to_dict()
    for i,r in enumerate(bam.fetch(until_eof=True)):
        if i >= 100000:
            break
        counts['records_sampled'] += 1
        counts['supplementary'] += r.is_supplementary
        counts['secondary'] += r.is_secondary
        counts['unmapped'] += r.is_unmapped
        counts['paired'] += r.is_paired
        counts['proper_pair'] += r.is_proper_pair
        counts['both_mates_mapped_different_contigs'] += bool(r.is_paired and not r.is_unmapped and not r.mate_is_unmapped and r.reference_id != r.next_reference_id)
        counts['with_intronic_N_CIGAR'] += any(op == 3 for op,n in (r.cigartuples or []))
        counts['with_soft_clip_CIGAR'] += any(op == 4 for op,n in (r.cigartuples or []))
        tags.update(tag for tag,value in r.get_tags())
full_qc = json.loads((out/'bam-qc.json').read_text())
rna_full = next(f for f in full_qc['files'] if f['kind'] == 'RNA')['flagstat']
disk = shutil.disk_usage(source)
res = {
    'generated_at': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'source': str(source/'RNA_TN26-279853.bam'),
    'pysam_version': pysam.__version__,
    'sample_scope': 'First 100,000 coordinate-ordered BAM records; limited tag/architecture audit, not a representative fusion discovery sample',
    'sample_counts': dict(counts),
    'sample_tag_counts': dict(tags),
    'full_file_flagstat_existing': rna_full,
    'header_programs': header.get('PG',[]),
    'header_sort_order': header.get('HD',{}),
    'resident_file_names': sorted(p.name for p in source.iterdir() if p.is_file()),
    'host_ram_bytes': int(subprocess.check_output(['sysctl','-n','hw.memsize']).strip()),
    'filesystem_available_bytes': disk.free,
    'references_downloaded': False,
    'alignment_or_fusion_calling_performed': False,
}
(out/'fusion-feasibility.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps({k:res[k] for k in ('sample_counts','sample_tag_counts','host_ram_bytes','filesystem_available_bytes')},indent=2))
