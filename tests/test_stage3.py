import struct
import unittest

from genopred_v2.stage3 import (
    _parse_gencode_hgnc_metadata_line,
    _parse_gtf_attributes,
    _unversion_ensembl_id,
)


class Stage3CompletionTests(unittest.TestCase):
    def test_bigwig_magic_constant_serializes_little_endian(self):
        self.assertEqual(bytes.fromhex("26fc8f88"), struct.pack("<I", 0x888FFC26))

    def test_gtf_attributes_are_parsed(self):
        attributes = _parse_gtf_attributes(
            'gene_id "ENSG000001.2"; gene_type "protein_coding"; gene_name "ABC";'
        )
        self.assertEqual("ENSG000001.2", attributes["gene_id"])
        self.assertEqual("ABC", attributes["gene_name"])

    def test_ensembl_unversioning_retains_par_y_suffix(self):
        self.assertEqual(
            "ENSG00000182378_PAR_Y",
            _unversion_ensembl_id("ENSG00000182378.14_PAR_Y"),
        )

    def test_hgnc_metadata_uses_identifier_column_not_symbol_column(self):
        self.assertEqual(
            ("ENST00000456328.2", "HGNC:37102"),
            _parse_gencode_hgnc_metadata_line(
                "ENST00000456328.2\tDDX11L1\tHGNC:37102\n"
            ),
        )


if __name__ == "__main__":
    unittest.main()
