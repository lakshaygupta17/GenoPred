import unittest
from pathlib import Path

from genopred_v2.provenance import load_source_manifest
from genopred_v2.schema import SchemaError, load_feature_schema, validate_feature_matrix


ROOT = Path(__file__).resolve().parents[1]


class ConfigurationTests(unittest.TestCase):
    def test_manifests_load_and_feature_sources_exist(self):
        manifest = load_source_manifest(ROOT / "configs" / "data_sources.json")
        features = load_feature_schema(ROOT / "configs" / "features.json")
        source_ids = {source.source_id for source in manifest.sources}
        referenced = {source_id for feature in features for source_id in feature.source_ids}
        self.assertEqual(set(), referenced - source_ids)

    def test_acquired_core_sources_are_recorded(self):
        manifest = load_source_manifest(ROOT / "configs" / "data_sources.json")
        statuses = {source.source_id: source.status for source in manifest.sources}
        for source_id in ("GSE149103", "GSE124229", "TCGA-PAAD"):
            self.assertEqual("downloaded_core", statuses[source_id])

    def test_feature_matrix_rejects_duplicate_gene_ids(self):
        with self.assertRaises(SchemaError):
            validate_feature_matrix(["KRAS", "KRAS"], [[1.0], [2.0]], expected_features=1)

    def test_feature_matrix_rejects_nan(self):
        with self.assertRaises(SchemaError):
            validate_feature_matrix(["KRAS"], [[float("nan")]], expected_features=1)


if __name__ == "__main__":
    unittest.main()
