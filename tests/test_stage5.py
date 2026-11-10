import gzip
import tempfile
import unittest
from pathlib import Path

import numpy as np

from genopred_v2.stage5 import (
    Stage5Error,
    _unique_liftover,
    aggregate_tcga_cnv,
    aggregate_tcga_mutations,
    collapse_patient_replicates,
    impute_with_missing_flag,
    modal_copy_baseline,
)


class FakeLiftOver:
    def __init__(self, mappings):
        self.mappings = mappings

    def convert_coordinate(self, chromosome, position):
        return self.mappings.get((chromosome, position), [])


class Stage5DatasetTests(unittest.TestCase):
    def test_tab_maf_deduplicates_cases_and_validates_schema(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'data/raw/TCGA-PAAD'
            source.mkdir(parents=True)
            path = source / 'sample.maf.gz'
            content = ('#version\nGene\tVariant_Classification\tTumor_Sample_Barcode\n'
                       'ENSG1.2\tMissense_Mutation\tTCGA-AB-1234-01A\n'
                       'ENSG1.2\tMissense_Mutation\tTCGA-AB-1234-01A\n'
                       'ENSG2.1\tSilent\tTCGA-AB-1234-01A\n')
            with gzip.open(path, 'wt') as f:
                f.write(content)
            cohort = [{'maf_file': path.name, 'case_id': 'TCGA-AB-1234'}]
            values, counts = aggregate_tcga_mutations(root, cohort, {'ENSG1':0,'ENSG2':1}, {'Missense_Mutation'})
            np.testing.assert_array_equal(values, [1, 0])
            np.testing.assert_array_equal(counts, [1, 1])
            with gzip.open(path, 'wt') as f:
                f.write(content.replace('\t', ','))
            with self.assertRaises(Stage5Error):
                aggregate_tcga_mutations(root, cohort, {'ENSG1':0}, {'Missense_Mutation'})

    def test_cnv_uses_case_baseline_and_excludes_sex_chromosomes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / 'data/raw/TCGA-PAAD'
            source.mkdir(parents=True)
            (source / 'cnv.tsv').write_text('gene_id\tcopy_number\nA\t4\nB\t4\nC\t8\nD\t2\nX\t1\n')
            genes = {g:i for i,g in enumerate('ABCDX')}
            identity = {g:{'chromosome':'chrX' if g=='X' else 'chr1'} for g in genes}
            amp, loss, missing, observed, qc = aggregate_tcga_cnv(root, [{'case_id':'case','cnv_file':'cnv.tsv'}], genes, identity)
            np.testing.assert_array_equal(amp[:4], [0,0,1,0])
            np.testing.assert_array_equal(loss[:4], [0,0,0,1])
            self.assertTrue(missing[4])
            self.assertEqual(0, observed[4])
            self.assertEqual(4, qc[0]['autosomal_modal_baseline'])
            self.assertEqual(4, modal_copy_baseline([2,4]))

    def test_patient_variability_excludes_within_patient_library_noise(self):
        profiles, patients = collapse_patient_replicates(np.array([[0,10,4,6]]), ['PT1_1','PT1_2','PT2_1','PT2_2'])
        np.testing.assert_array_equal(profiles, [[5,5]])
        self.assertEqual(['PT1','PT2'], patients)

    def test_liftover_rejects_expanded_and_reversed_intervals(self):
        for end in [100000, 10]:
            converter = FakeLiftOver({('chr1',0):[('chr1',20,'+',1)], ('chr1',9):[('chr1',end,'+',1)]})
            self.assertIsNone(_unique_liftover(converter, 'chr1',0,10))

    def test_imputation_adds_explicit_missing_flag(self):
        values, flags, median = impute_with_missing_flag(
            np.asarray([1.0, np.nan, 3.0]), np.asarray([False, False, False])
        )
        np.testing.assert_array_equal(values, np.asarray([1.0, 2.0, 3.0]))
        np.testing.assert_array_equal(flags, np.asarray([0.0, 1.0, 0.0]))
        self.assertEqual(2.0, median)

    def test_liftover_requires_unique_consistent_endpoints(self):
        converter = FakeLiftOver(
            {
                ("chr1", 10): [("chr1", 20.0, "+", 1.0)],
                ("chr1", 19): [("chr1", 29.0, "+", 1.0)],
            }
        )
        self.assertEqual(("chr1", 20, 30), _unique_liftover(converter, "chr1", 10, 20))

    def test_liftover_rejects_cross_chromosome_mapping(self):
        converter = FakeLiftOver(
            {
                ("chr1", 10): [("chr1", 20.0, "+", 1.0)],
                ("chr1", 19): [("chr2", 29.0, "+", 1.0)],
            }
        )
        self.assertIsNone(_unique_liftover(converter, "chr1", 10, 20))


if __name__ == "__main__":
    unittest.main()
