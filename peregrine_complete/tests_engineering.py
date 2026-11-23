"""Offline contract tests for engineering entry point."""
import json
import tempfile
import unittest
from pathlib import Path
from peregrine_engineering import validate_links, validate_manifest, main

class EngineeringTests(unittest.TestCase):
    def test_valid_and_empty_links(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'links'; p.write_text('enh\tPANTHER1\t55\t1\n')
            self.assertEqual(validate_links(p,'1')['rows'],1)
            p.write_text(''); self.assertEqual(validate_links(p)['rows'],0)
    def test_invalid_links(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td)/'links'; p.write_text('enh\tgene\t55\n')
            with self.assertRaises(ValueError): validate_links(p)
            p.write_text('enh\tgene\t55\t2\n')
            with self.assertRaises(ValueError): validate_links(p,'1')
    def test_manifest(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td); (root/'s.mango').write_text('chr1\t1\t2\tchr1\t3\t4\t1\t0.01\n')
            p=root/'samples.tsv'; p.write_text('sample_id\tcell_type\ttarget_protein\tmango_file\nS1\tK562\tPOLR2A\ts.mango\n')
            self.assertEqual(len(validate_manifest(p)),1)
            p.write_text(p.read_text()+'S1\tK562\tPOLR2A\ts.mango\n')
            with self.assertRaises(ValueError): validate_manifest(p)
    def test_failed_run_writes_report(self):
        with tempfile.TemporaryDirectory() as td:
            out=Path(td)/'results'
            self.assertEqual(main(['validate','--output',str(out)]),1)
            self.assertEqual(json.loads((out/'engineering/report.json').read_text())['status'],'failed')

if __name__=='__main__': unittest.main()
