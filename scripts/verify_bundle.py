#!/usr/bin/env python3
"""Verify exact repository file bytes against the supplied integrity manifest."""
import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    sha = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            sha.update(chunk)
    return sha.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = root / 'manifests/repository-files.sha256.json'
    rows = json.loads(manifest.read_text())['files']
    failures = []
    for row in rows:
        relative = Path(row['path'])
        if relative.is_absolute() or '..' in relative.parts:
            failures.append((row['path'], 'unsafe manifest path'))
            continue
        path = root / relative
        if not path.is_file():
            failures.append((row['path'], 'missing'))
        elif path.stat().st_size != row['size_bytes']:
            failures.append((row['path'], 'size mismatch'))
        elif digest(path) != row['sha256']:
            failures.append((row['path'], 'SHA256 mismatch'))
    for name, reason in failures:
        print(f'FAIL {reason}: {name}')
    print(f'{len(rows) - len(failures)}/{len(rows)} manifest-listed repository files verified.')
    print('The manifest excludes itself and .git/. Extra unlisted files are not checked.')
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
