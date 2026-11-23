"""Offline smoke and contract tests, runnable with python -m unittest tests_preprocessing."""
import gzip
import tempfile
import unittest
from pathlib import Path
from preprocess.auto import gtf_to_tables, prepare_annotations, prepare_eqtl, prepare_tad

class TestWorkflow(unittest.TestCase):
    def test_offline_components(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / 'data' / 'raw' / 'gtex' / 'GTEx_Analysis_v7_eQTL').mkdir(parents=True)
            (root / 'data' / 'raw' / 'encode').mkdir()
            gtf = root / 'sample.gtf.gz'
            with gzip.open(gtf, 'wt') as out:
                for strand, gene, start, end in [('+', 'ENSG000001', 1000, 2000), ('-', 'ENSG000002', 3000, 4000)]:
                    attrs = f'gene_id "{gene}"; gene_biotype "protein_coding";'
                    out.write(f'1\tENSEMBL\tgene\t{start}\t{end}\t.\t{strand}\t.\t{attrs}\n')
                    out.write(f'1\tENSEMBL\texon\t{start}\t{start+99}\t.\t{strand}\t.\t{attrs}\n')
            outputs = prepare_annotations(root, gtf=gtf, download=False)
            self.assertIn('chr1\t4000\t4600\tENSG000002', outputs['promoters'].read_text())
            self.assertIn('chr1\t400\t1000\tENSG000001', outputs['promoters'].read_text())
            self.assertEqual(len(outputs['exons'].read_text().splitlines()), 2)
            self.assertEqual(len(outputs['tss'].read_text().splitlines()), 2)
            eqtl = root / 'data' / 'raw' / 'gtex' / 'GTEx_Analysis_v7_eQTL' / 'Liver.signif_variant_gene_pairs.txt.gz'
            with gzip.open(eqtl, 'wt') as out:
                out.write('variant_id\tgene_id\tpval_nominal\n1_123_A_G_b37\tENSG000001\t0.001\n')
            self.assertEqual(prepare_eqtl(root), 1)
            self.assertIn('chr1\t123\t123', next((root / 'data' / 'processed' / 'eqtl').glob('*.bed')).read_text())
            (root / 'data' / 'raw' / 'encode' / 'ENCFF274VJU.bed').write_text('1\t100\t200\n1\t-1\t200\n')
            self.assertEqual(prepare_tad(root), 1)
            self.assertEqual((root / 'data' / 'processed' / 'tad' / 'ENCFF274VJU.bed').read_text(), 'chr1\t100\t200\n')

if __name__ == '__main__':
    unittest.main()
