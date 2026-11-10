import unittest

from genopred_v2.acquisition import (
    _gdc_selection,
    _geo_selection,
    _parse_soft_records,
    _values_with_prefix,
)


class AcquisitionTests(unittest.TestCase):
    def test_soft_parser_preserves_repeated_fields(self):
        text = """^SAMPLE = GSM1
!Sample_title = PANC-1 ATAC-seq Rep1
!Sample_supplementary_file_1 = ftp://example/one.bw
!Sample_supplementary_file_2 = ftp://example/two.bw
"""
        records = _parse_soft_records(text)
        self.assertEqual(1, len(records))
        self.assertEqual("SAMPLE", records[0][0])
        self.assertEqual("GSM1", records[0][1])
        files = _values_with_prefix(records[0][2], "Sample_supplementary_file")
        self.assertEqual(["ftp://example/one.bw", "ftp://example/two.bw"], files)

    def test_panc1_atac_is_core(self):
        selection, _ = _geo_selection(
            "GSE149103", "GSM4490514", "PANC-1 ATAC-seq Rep1", "PANC-1_ATAC.bw"
        )
        self.assertEqual("core", selection)

    def test_capan1_is_validation(self):
        selection, _ = _geo_selection(
            "GSE149103", "GSM4490536", "Capan-1 ATAC-seq Rep1", "Capan-1_ATAC.bw"
        )
        self.assertEqual("validation", selection)

    def test_patient_atac_count_matrix_is_core(self):
        selection, _ = _geo_selection(
            "GSE124229", "GSE124229", "PDAC ATAC", "GSE124229_Atlas_counts.txt.gz"
        )
        self.assertEqual("core", selection)

    def test_patient_bigwig_is_optional_validation(self):
        selection, _ = _geo_selection(
            "GSE124229", "GSM3526038", "PT12_1", "GSM3526038_sample.bw"
        )
        self.assertEqual("validation", selection)

    def test_only_primary_cnv_workflow_is_core(self):
        self.assertEqual("core", _gdc_selection("Copy Number Variation", "Gene Level Copy Number", "ASCAT3")[0])
        self.assertEqual("validation", _gdc_selection("Copy Number Variation", "Gene Level Copy Number", "ABSOLUTE LiftOver")[0])
        self.assertEqual("supporting", _gdc_selection("Copy Number Variation", "Gene Level Copy Number", "ASCAT2")[0])


if __name__ == "__main__":
    unittest.main()
