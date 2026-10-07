#!/usr/bin/env python3
"""Check published aggregate arithmetic and local Markdown destinations.

Uses only Python's standard library. No sequences, reads, network requests or
patient-level processing. Passing does not independently reproduce read counts.
"""
import argparse
import json
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit


class ReviewError(ValueError):
    pass


def population_checks(p):
    keys = ['selected_sites', 'original_joint_callable', 'original_dna_dominant',
            'original_same90', 'original_same98', 'new_dominant_callable',
            'new_same90_fixed', 'new_same98_fixed', 'new_depth_lost',
            'new_callable_below90', 'strong_discrepancies_before', 'strong_discrepancies_after']
    if not isinstance(p, dict):
        raise ReviewError('population_panel must be an object')
    for key in keys:
        if key not in p or type(p[key]) is not int or p[key] < 0:
            raise ReviewError(f'{key} must be a nonnegative integer')
    n = p['original_dna_dominant']
    if n == 0:
        raise ReviewError('fixed denominator must be positive')
    checks = {
      'selection_and_original_depth': p['selected_sites'] >= p['original_joint_callable'] >= n,
      'original_threshold_nesting': n >= p['original_same90'] >= p['original_same98'],
      'new_threshold_nesting': p['new_dominant_callable'] >= p['new_same90_fixed'] >= p['new_same98_fixed'],
      'callable_plus_dropout_equals_fixed_denominator': p['new_dominant_callable'] + p['new_depth_lost'] == n,
      'new_callable_partition': p['new_same90_fixed'] + p['new_callable_below90'] == p['new_dominant_callable'],
      'discrepancy_counts_are_bounded': p['strong_discrepancies_before'] <= n and p['strong_discrepancies_after'] <= n,
    }
    t = p.get('transitions')
    allowed = {'retain90__retain90', 'retain90__below20', 'below90__below90',
               'below90__retain90', 'retain90__below90'}
    if not isinstance(t, dict) or set(t) != allowed or any(type(x) is not int or x < 0 for x in t.values()):
        raise ReviewError('transition table is missing, malformed or has unrecognized categories')
    checks.update({
      'transition_total_is_fixed_denominator': sum(t.values()) == n,
      'transition_original_retention': t['retain90__retain90'] + t['retain90__below20'] + t['retain90__below90'] == p['original_same90'],
      'transition_new_retention': t['retain90__retain90'] + t['below90__retain90'] == p['new_same90_fixed'],
      'transition_dropout': t['retain90__below20'] == p['new_depth_lost'],
      'transition_remaining_below90': t['below90__below90'] + t['retain90__below90'] == p['new_callable_below90'],
    })
    failures = [name for name, ok in checks.items() if not ok]
    if failures:
        raise ReviewError('inconsistent aggregates: ' + ', '.join(failures))
    return {'checks_passed': len(checks), 'fixed_denominator': n,
            'original_same90_percent': 100*p['original_same90']/n,
            'original_same98_percent': 100*p['original_same98']/n,
            'new_same90_fixed_percent': 100*p['new_same90_fixed']/n,
            'new_same98_fixed_percent': 100*p['new_same98_fixed']/n}


def local_links(root):
    root = Path(root).resolve()
    failures, checked = [], 0
    # Curated release documents use inline Markdown destinations and reference
    # definitions. External URLs and within-document anchors are not fetched.
    inline = re.compile(r'\]\(\s*(<[^>]+>|[^\s)]+)(?:\s+["\'][^\n]*?["\'])?\s*\)')
    reference = re.compile(r'^\s*\[[^\]]+\]:\s*(<[^>]+>|\S+)', re.M)
    for doc in root.rglob('*.md'):
        if any(part in {'.git', '__pycache__', '.pytest_cache'} for part in doc.relative_to(root).parts):
            continue
        text = doc.read_text(encoding='utf-8')
        text = re.sub(r'```.*?```', '', text, flags=re.S)
        for match in list(inline.finditer(text)) + list(reference.finditer(text)):
            value = match.group(1).strip('<>')
            split = urlsplit(value)
            if split.scheme or split.netloc or not split.path:
                continue
            checked += 1
            relative = unquote(split.path)
            if relative.startswith('/') or '\\' in relative:
                failures.append(f'{doc.relative_to(root)}: nonportable link {value}')
                continue
            target = (doc.parent / relative).resolve()
            if not target.is_relative_to(root):
                failures.append(f'{doc.relative_to(root)}: link escapes release {value}')
            elif not target.exists():
                failures.append(f'{doc.relative_to(root)}: missing destination {value}')
    return {'checked_local_destinations': checked, 'failures': failures,
            'limits': 'Checks file/directory existence, not heading anchors, HTML links or remote availability.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    try:
        data = json.loads((args.root/'data/review-summary.json').read_text(encoding='utf-8'))
        result = {'aggregates': population_checks(data['population_panel']), 'links': local_links(args.root),
                  'scope': 'Published-summary consistency and local navigation only; no independent patient-read verification.'}
    except (OSError, UnicodeError, json.JSONDecodeError, KeyError, ReviewError) as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2))
    return 1 if result['links']['failures'] else 0


if __name__ == '__main__':
    raise SystemExit(main())
