"""Synthetic, non-patient regression and negative cases for public helpers."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from check_release import ReviewError, local_links, population_checks
from synthetic_counting import parse_line
from verify_bundle import VerificationError, verify
from verify_inputs import verify_inputs


PUBLISHED = {
    'selected_sites': 165782, 'original_joint_callable': 4015,
    'original_dna_dominant': 2873, 'original_same90': 2860,
    'original_same98': 2760, 'new_dominant_callable': 2814,
    'new_same90_fixed': 2805, 'new_same98_fixed': 2721,
    'new_depth_lost': 59, 'new_callable_below90': 9,
    'strong_discrepancies_before': 6, 'strong_discrepancies_after': 5,
    'transitions': {'retain90__retain90': 2800, 'retain90__below20': 59,
                    'below90__below90': 8, 'below90__retain90': 5,
                    'retain90__below90': 1},
}


class AggregateTests(unittest.TestCase):
    def test_published_count_arithmetic(self):
        r = population_checks(PUBLISHED)
        self.assertEqual(r['fixed_denominator'], 2873)
        self.assertAlmostEqual(r['original_same90_percent'], 99.5475, places=3)
        self.assertAlmostEqual(r['new_same90_fixed_percent'], 97.6331, places=3)

    def test_fraction_only_count_cannot_replace_depth_qualified_count(self):
        p = copy.deepcopy(PUBLISHED)
        p['new_same90_fixed'] = 2860  # Includes shallow sites in the historical pitfall.
        with self.assertRaises(ReviewError):
            population_checks(p)

    def test_dropping_coverage_failures_from_denominator_is_detected(self):
        p = copy.deepcopy(PUBLISHED)
        p['original_dna_dominant'] = 2814
        with self.assertRaises(ReviewError):
            population_checks(p)

    def test_transition_with_unchanged_total_but_wrong_destination_is_detected(self):
        p = copy.deepcopy(PUBLISHED)
        p['transitions']['retain90__below20'] -= 1
        p['transitions']['retain90__below90'] += 1
        with self.assertRaises(ReviewError):
            population_checks(p)

    def test_zero_denominator_and_boolean_counts_are_rejected(self):
        for key, value in [('original_dna_dominant', 0), ('new_depth_lost', True)]:
            p = copy.deepcopy(PUBLISHED)
            p[key] = value
            with self.assertRaises(ReviewError):
                population_checks(p)


class ManifestTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name) / 'release'
        self.root.mkdir()
        (self.root/'manifests').mkdir()
        (self.root/'example.txt').write_bytes(b'abc\n')
        self.row = {'path': 'example.txt', 'size_bytes': 4,
                    'sha256': hashlib.sha256(b'abc\n').hexdigest()}
        self.manifest = self.root/'manifests/repository-files.sha256.json'
        self.write_rows([self.row])

    def tearDown(self):
        self.tmp.cleanup()

    def write_rows(self, rows):
        self.manifest.write_text(json.dumps({'files': rows}))

    def test_clean_and_ignored_python_cache(self):
        (self.root/'__pycache__').mkdir()
        (self.root/'__pycache__/ignored.pyc').write_bytes(b'cache')
        self.assertEqual(verify(self.root)['failures'], [])

    def test_same_size_corruption_is_detected(self):
        (self.root/'example.txt').write_bytes(b'xyz\n')
        self.assertTrue(any('SHA256 mismatch' in x for x in verify(self.root)['failures']))

    def test_extra_file_is_not_silently_ignored(self):
        (self.root/'unintended.txt').write_text('synthetic test')
        self.assertTrue(any('unlisted file' in x for x in verify(self.root)['failures']))

    def test_missing_file_is_detected(self):
        (self.root/'example.txt').unlink()
        self.assertTrue(any('missing' in x for x in verify(self.root)['failures']))

    def test_empty_duplicate_and_invalid_metadata(self):
        for rows in [[], [self.row, self.row], [{**self.row, 'size_bytes': True}],
                     [{**self.row, 'sha256': 'not-a-hash'}]]:
            self.write_rows(rows)
            with self.assertRaises(VerificationError):
                verify(self.root)

    def test_traversal_and_absolute_paths_are_rejected(self):
        for path in ['../example.txt', '/example.txt', 'a/../example.txt', 'a\\example.txt']:
            self.write_rows([{**self.row, 'path': path}])
            with self.assertRaises(VerificationError):
                verify(self.root)

    def test_symlink_file_and_symlink_parent_are_rejected(self):
        outside = Path(self.tmp.name)/'outside'
        outside.mkdir()
        (outside/'data.txt').write_bytes(b'abc\n')
        (self.root/'link').symlink_to(outside, target_is_directory=True)
        self.write_rows([{**self.row, 'path': 'link/data.txt'}])
        with self.assertRaises(VerificationError):
            verify(self.root)


class NavigationTests(unittest.TestCase):
    def test_relative_links_spaces_and_external_links(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root/'a b.txt').write_text('synthetic')
            (root/'README.md').write_text('[ok](a%20b.txt) [remote](https://example.org/no-fetch)')
            self.assertEqual(local_links(root)['failures'], [])
            (root/'README.md').write_text('[missing](missing.txt) [escape](../outside.txt)')
            self.assertEqual(len(local_links(root)['failures']), 2)


class InputManifestTests(unittest.TestCase):
    def test_valid_and_same_size_corruption(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            (root/'input.txt').write_bytes(b'abc\n')
            row = {'name':'input.txt','size_bytes':4,
                   'sha256':hashlib.sha256(b'abc\n').hexdigest(), 'included_in_repository':True}
            manifest = root/'inputs.json'
            manifest.write_text(json.dumps({'files':[row], 'total_files':1,'total_bytes':4}))
            self.assertEqual(verify_inputs(root, manifest, True)['failures'], [])
            (root/'input.txt').write_bytes(b'xyz\n')
            self.assertTrue(verify_inputs(root, manifest)['failures'])

    def test_external_absolute_path_and_duplicate_names_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            row = {'name':'input.txt','size_bytes':0,'sha256':hashlib.sha256(b'').hexdigest()}
            manifest = root/'inputs.json'
            for rows in [[{**row,'name':str(root/'outside.txt')}],[row,row]]:
                manifest.write_text(json.dumps({'files':rows}))
                with self.assertRaises(VerificationError):
                    verify_inputs(root, manifest)

    def test_empty_selection_and_wrong_manifest_total_are_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root = Path(d)
            row = {'name':'input.txt','size_bytes':0,'sha256':hashlib.sha256(b'').hexdigest(),'included_in_repository':False}
            manifest = root/'inputs.json'
            manifest.write_text(json.dumps({'files':[row]}))
            with self.assertRaises(VerificationError):
                verify_inputs(root, manifest, True)
            manifest.write_text(json.dumps({'files':[row],'total_bytes':1}))
            with self.assertRaises(VerificationError):
                verify_inputs(root, manifest)


class SyntheticCountingTests(unittest.TestCase):
    marker = {'chrom': 'synthetic', 'pos1': 101, 'ref': 'A', 'alt': 'G'}

    def line(self, bases, names):
        n = len(bases)
        return '\t'.join(['synthetic', '101', 'N', str(n), bases, 'D'*n, ']'*n,
                          ','.join(names), ','.join(['31']*n)]) + '\n'

    def test_overlapping_mates_count_once_and_conflict_is_excluded(self):
        result = parse_line(self.line('AaGgAg', ['ref','ref','alt','alt','conflict','conflict']), self.marker)
        self.assertEqual((result['read_ref'], result['read_alt']), (3, 3))
        self.assertEqual((result['fragment_ref'], result['fragment_alt'], result['fragment_discordant']), (1, 1, 1))
        self.assertEqual(result['fragment_depth_refalt'], 2)

    def test_deletion_and_refskip_are_not_base_coverage(self):
        result = parse_line(self.line('*<>', ['deletion','skip-forward','skip-reverse']), self.marker)
        self.assertEqual(result['fragment_nonACGT'], 3)
        self.assertEqual(result['fragment_depth_refalt'], 0)
        self.assertIsNone(result['fragment_ALT_fraction'])

    def test_zero_depth_and_malformed_parallel_fields(self):
        result = parse_line('synthetic\t101\tN\t0\t*\t*\t*\t*\t*\n', self.marker)
        self.assertEqual(result['fragment_depth_refalt'], 0)
        with self.assertRaises(ValueError):
            parse_line('synthetic\t101\tN\t2\tAG\tDD\t]]\tonly-one-name\t31,31\n', self.marker)


if __name__ == '__main__':
    unittest.main()
