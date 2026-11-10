import gzip
import tempfile
import unittest
from pathlib import Path

from genopred_v2.qc import (
    classify_tcga_file,
    infer_gse149103_sample,
    read_header,
    schema_signature,
    _sample_associations,
)


class Stage3QCTests(unittest.TestCase):
    def test_tcga_modalities_are_classified(self):
        self.assertEqual(
            "rna_expression",
            classify_tcga_file("x.rna_seq.augmented_star_gene_counts.tsv"),
        )
        self.assertEqual(
            "gene_level_cnv",
            classify_tcga_file("TCGA-PAAD.x.ascat3.gene_level_copy_number.v36.tsv"),
        )
        self.assertEqual("somatic_mutation", classify_tcga_file("x.maf.gz"))

    def test_gse149103_filename_metadata(self):
        accession, cell_line, assay = infer_gse149103_sample(
            "GSM4490518_PANC-1_1_H3K27ac.bw"
        )
        self.assertEqual(("GSM4490518", "PANC-1", "H3K27ac ChIP-seq"), (accession, cell_line, assay))

    def test_gzip_tabular_header_is_read_after_comments(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "table.tsv.gz"
            with gzip.open(path, "wt", encoding="utf-8") as handle:
                handle.write("# build: GRCh38\na\tb\n1\t2\n")
            comments, columns, first_row = read_header(path)
        self.assertEqual(["# build: GRCh38"], comments)
        self.assertEqual(["a", "b"], columns)
        self.assertEqual(["1", "2"], first_row)

    def test_schema_signature_is_stable_and_order_sensitive(self):
        self.assertEqual(schema_signature(["a", "b"]), schema_signature(["a", "b"]))
        self.assertNotEqual(schema_signature(["a", "b"]), schema_signature(["b", "a"]))

    def test_exact_gdc_entity_is_joined_to_sample_metadata(self):
        hit = {
            "cases": [{"submitter_id": "TCGA-AB-1234", "samples": [{
                "submitter_id": "TCGA-AB-1234-01A",
                "sample_type": "Primary Tumor",
                "tissue_type": "Tumor",
            }]}],
            "associated_entities": [{
                "entity_submitter_id": "TCGA-AB-1234-01A-11R-0000-00",
                "entity_type": "aliquot",
            }],
        }
        association = _sample_associations(hit)[0]
        self.assertEqual("TCGA-AB-1234", association["subject_id"])
        self.assertEqual("Primary Tumor", association["sample_type"])


if __name__ == "__main__":
    unittest.main()
