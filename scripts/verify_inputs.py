#!/usr/bin/env python3
"""Verify manifest-listed original bytes without downloading or modifying.

Input contents are never printed. Other data-directory files are not examined.
Passing verifies bytes against a manifest, not scientific validity.
"""
import argparse
import json
from pathlib import Path, PurePosixPath
import re
import sys

from verify_bundle import VerificationError, digest, safe_path


def verify_inputs(data_dir, manifest, small_only=False):
    root = Path(data_dir).resolve()
    if not root.is_dir():
        raise VerificationError('data directory does not exist')
    try:
        doc = json.loads(Path(manifest).read_text(encoding='utf-8'))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise VerificationError(f'cannot read input manifest: {exc}') from exc
    records = doc.get('files') if isinstance(doc, dict) else None
    if not isinstance(records, list) or not records:
        raise VerificationError('input manifest files must be a nonempty list')
    seen = set()
    for row in records:
        if not isinstance(row, dict) or not {'name','size_bytes','sha256'} <= row.keys():
            raise VerificationError('malformed input manifest record')
        name = row['name']
        safe_path(root, name)
        if len(PurePosixPath(name).parts) != 1 or name in seen:
            raise VerificationError('input names must be unique basenames')
        seen.add(name)
        if type(row['size_bytes']) is not int or row['size_bytes'] < 0:
            raise VerificationError(f'invalid input size: {name!r}')
        if not isinstance(row['sha256'], str) or not re.fullmatch(r'[0-9a-f]{64}', row['sha256']):
            raise VerificationError(f'invalid input SHA256: {name!r}')
        if small_only and type(row.get('included_in_repository')) is not bool:
            raise VerificationError('small-only selection needs explicit boolean included_in_repository values')
    if 'total_files' in doc and (type(doc['total_files']) is not int or doc['total_files'] != len(records)):
        raise VerificationError('manifest total_files disagrees with file list')
    if 'total_bytes' in doc and (type(doc['total_bytes']) is not int or doc['total_bytes'] != sum(x['size_bytes'] for x in records)):
        raise VerificationError('manifest total_bytes disagrees with file list')
    if small_only:
        records = [x for x in records if x['included_in_repository']]
    if not records:
        raise VerificationError('selected input set is empty')
    failures = []
    for row in records:
        path = safe_path(root, row['name'])
        try:
            if not path.is_file():
                failures.append(f'missing/not regular: {row["name"]}')
            elif path.stat().st_size != row['size_bytes']:
                failures.append(f'size mismatch: {row["name"]}')
            elif digest(path) != row['sha256']:
                failures.append(f'SHA256 mismatch: {row["name"]}')
        except OSError as exc:
            failures.append(f'unreadable: {row["name"]}: {exc}')
    return {'inputs_checked':len(records), 'failures':failures,
            'scope':'Manifest-listed file integrity only; unlisted files and scientific validity are not checked.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, default=Path(__file__).resolve().parents[1] / 'manifests/raw-inputs.json')
    parser.add_argument('--small-only', action='store_true', help='Verify only the three small inputs included in the repository.')
    args = parser.parse_args()
    try:
        result = verify_inputs(args.data_dir, args.manifest, args.small_only)
    except VerificationError as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2))
    return 1 if result['failures'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
