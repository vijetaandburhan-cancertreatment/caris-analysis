#!/usr/bin/env python3
"""Hash the nine original Caris inputs. Does not download or modify anything."""
import argparse
import hashlib
import json
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir', type=Path, required=True)
    parser.add_argument('--manifest', type=Path, default=Path(__file__).resolve().parents[1] / 'manifests/raw-inputs.json')
    parser.add_argument('--small-only', action='store_true', help='Verify only the three small inputs included in the repository.')
    args = parser.parse_args()
    records = json.loads(args.manifest.read_text())['files']
    if args.small_only:
        records = [r for r in records if r['included_in_repository']]
    failed = 0
    for record in records:
        path = args.data_dir / record['name']
        if not path.is_file():
            print('MISSING', record['name'], flush=True)
            failed += 1
            continue
        if path.stat().st_size != record['size_bytes']:
            print('SIZE MISMATCH', record['name'], flush=True)
            failed += 1
            continue
        valid = sha256(path) == record['sha256']
        print('OK' if valid else 'HASH MISMATCH', record['name'], flush=True)
        failed += not valid
    print(f'{len(records) - failed}/{len(records)} inputs verified.')
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
