import unittest

from genopred_v2.stage4 import classify_label, resolve_gene_symbol


class Stage4LabelTests(unittest.TestCase):
    def test_two_external_families_qualify_positive(self):
        self.assertEqual("positive", classify_label(["CPTAC", "ICGC_PACA_AU"], 2))

    def test_duplicate_support_within_one_family_does_not_qualify(self):
        self.assertEqual("unlabelled", classify_label(["ICGC_PACA_AU", "ICGC_PACA_AU"], 2))

    def test_tcga_only_evidence_does_not_qualify(self):
        self.assertEqual("unlabelled", classify_label([], 2))

    def test_symbol_resolution_rejects_ambiguous_symbols(self):
        status, row = resolve_gene_symbol(
            "DUP", {"DUP": [{"ensembl_gene_id": "ENSG1"}, {"ensembl_gene_id": "ENSG2"}]}
        )
        self.assertEqual("ambiguous", status)
        self.assertIsNone(row)

    def test_symbol_resolution_accepts_unique_symbol(self):
        expected = {"ensembl_gene_id": "ENSG1"}
        status, row = resolve_gene_symbol("ABC", {"ABC": [expected]})
        self.assertEqual("unique", status)
        self.assertEqual(expected, row)


if __name__ == "__main__":
    unittest.main()
