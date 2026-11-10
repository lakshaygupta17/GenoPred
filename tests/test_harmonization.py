import unittest

from genopred_v2.harmonization import (
    GenomeBuildMismatchError,
    require_compatible_build,
    requires_harmonization,
)


class HarmonizationTests(unittest.TestCase):
    def test_hg19_is_recognized_as_grch37(self):
        self.assertFalse(requires_harmonization("hg19", "GRCh37"))

    def test_hg19_cannot_silently_mix_with_grch38(self):
        with self.assertRaises(GenomeBuildMismatchError):
            require_compatible_build("hg19", "GRCh38", "GSE149103")

    def test_matching_grch38_is_allowed(self):
        require_compatible_build("GRCh38", "GRCh38", "TCGA-PAAD")

    def test_unknown_build_is_rejected(self):
        with self.assertRaises(GenomeBuildMismatchError):
            require_compatible_build("unknown", "GRCh38", "GSE185069")


if __name__ == "__main__":
    unittest.main()

