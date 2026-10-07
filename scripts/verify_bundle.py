#!/usr/bin/env python3
"""Verify a release manifest and detect unlisted files; never downloads data.

Checks byte integrity, not provenance or scientific correctness. The manifest
itself and Git/Python cache directories are intentionally excluded.
"""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import sys

IGNORED_DIRECTORIES = {'.git', '__pycache__', '.pytest_cache'}


class VerificationError(ValueError):
    pass


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def safe_path(root, name):
    if not isinstance(name, str) or not name or '\\' in name or '\x00' in name or ':' in name:
        raise VerificationError('invalid relative path')
    rel = PurePosixPath(name)
    if rel.is_absolute() or '..' in rel.parts or name != rel.as_posix() or rel == PurePosixPath('.'):
        raise VerificationError(f'unsafe or noncanonical manifest path: {name!r}')
    if any(part in IGNORED_DIRECTORIES for part in rel.parts):
        raise VerificationError(f'excluded directory appears in manifest: {name!r}')
    path = root.joinpath(*rel.parts)
    current = root
    for part in rel.parts:
        current = current / part
        if current.is_symlink():
            raise VerificationError(f'symlink not allowed: {name!r}')
    if not path.resolve().is_relative_to(root):
        raise VerificationError(f'path escapes release root: {name!r}')
    return path


def verify(root, manifest_name='manifests/repository-files.sha256.json'):
    root = Path(root).resolve()
    if not root.is_dir():
        raise VerificationError('release root is not a directory')
    manifest = safe_path(root, manifest_name)
    try:
        document = json.loads(manifest.read_text(encoding='utf-8'))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise VerificationError(f'cannot read manifest: {exc}') from exc
    rows = document.get('files') if isinstance(document, dict) else None
    if not isinstance(rows, list) or not rows:
        raise VerificationError('manifest files must be a nonempty list')
    seen, failures = set(), []
    for row in rows:
        if not isinstance(row, dict) or not {'path', 'size_bytes', 'sha256'} <= row.keys():
            raise VerificationError('malformed manifest record')
        name, size, expected = row['path'], row['size_bytes'], row['sha256']
        path = safe_path(root, name)
        if name == manifest_name or name in seen:
            raise VerificationError(f'duplicate or self-referential manifest entry: {name!r}')
        seen.add(name)
        if type(size) is not int or size < 0 or not isinstance(expected, str) or not re.fullmatch(r'[0-9a-f]{64}', expected):
            raise VerificationError(f'invalid size/hash metadata: {name!r}')
        try:
            if not path.is_file():
                failures.append(f'missing/not regular: {name}')
            elif path.stat().st_size != size:
                failures.append(f'size mismatch: {name}')
            elif digest(path) != expected:
                failures.append(f'SHA256 mismatch: {name}')
        except OSError as exc:
            failures.append(f'unreadable: {name}: {exc}')
    actual = set()
    for path in root.rglob('*'):
        rel = path.relative_to(root)
        if any(part in IGNORED_DIRECTORIES for part in rel.parts):
            continue
        if path.is_symlink():
            failures.append(f'symlink not allowed: {rel.as_posix()}')
        elif path.is_file() and rel.as_posix() != manifest_name:
            actual.add(rel.as_posix())
    failures.extend(f'unlisted file: {name}' for name in sorted(actual - seen))
    return {'manifest_entries': len(rows), 'failures': failures,
            'scope': 'Byte-integrity and inventory check only; not provenance or scientific validation.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--manifest', default='manifests/repository-files.sha256.json')
    args = parser.parse_args()
    try:
        result = verify(args.root, args.manifest)
    except VerificationError as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2))
    return 1 if result['failures'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
