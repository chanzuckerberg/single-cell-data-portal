import unittest

from backend.wmg.data.tissue_mapper import TissueMapper
from backend.wmg.pipeline.primary_filter_dimensions import order_tissues


class TissueMapperTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tissue_mapper = TissueMapper()

    def test__high_level_tissue_retrieval_exists(self):
        low_level_tissue = "UBERON:0008951"  # lef lung lobe
        expected_high_level_tissue = "UBERON:0002048"  # lung
        self.assertEqual(self.tissue_mapper.get_high_level_tissue(low_level_tissue), expected_high_level_tissue)

    def test__high_level_tissue_retrieval_does_not_exist(self):
        low_level_tissue = "UBERON:noId"
        expected_high_level_tissue = "UBERON:noId"
        self.assertEqual(self.tissue_mapper.get_high_level_tissue(low_level_tissue), expected_high_level_tissue)

    def test__high_level_tissue_retrieval_suffix(self):
        low_level_tissue = "UBERON:0008951 (organoid)"  # lef lung lobe
        expected_high_level_tissue = "UBERON:0008951 (organoid)"  # lung
        self.assertEqual(self.tissue_mapper.get_high_level_tissue(low_level_tissue), expected_high_level_tissue)

    def test__making_ontology_id_writable(self):
        tissue = "UBERON_0008951"
        expected_tissue = "UBERON:0008951"

        self.assertEqual(self.tissue_mapper.reformat_ontology_term_id(tissue, to_writable=True), expected_tissue)

    def test__making_ontology_id_readable(self):
        tissue = "UBERON:0008951"
        expected_tissue = "UBERON_0008951"

        self.assertEqual(self.tissue_mapper.reformat_ontology_term_id(tissue, to_writable=False), expected_tissue)

    def test__get_label_from_id(self):
        tissue = "UBERON:0008951"
        expected_label = "left lung lobe"

        self.assertEqual(self.tissue_mapper.get_label_from_writable_id(tissue), expected_label)


class TissueOrderTest(unittest.TestCase):

    """
    Tests function used for ordering final tissue list sent to frontend
    """

    def test__tissue_ordering(self):
        unordered_tissues = [
            "UBERON:0002048 (organoid)",  # lung (organoid)
            "UBERON:0002106",  # spleen
            "UBERON:0002371",  # bone marrow
            "UBERON:0002107",  # liver
            "UBERON:0002113",  # kidney
            "UBERON:0000178",  # blood
            "UBERON:0002240",  # spinal cord
            "UBERON:0002048",  # lung
            "UBERON:0000955",  # brain
        ]

        expected_ordered_tissues = [
            "UBERON:0000178",  # blood
            "UBERON:0002048",  # lung
            "UBERON:0002106",  # spleen
            "UBERON:0002371",  # bone marrow
            "UBERON:0002107",  # liver
            "UBERON:0002113",  # kidney
            "UBERON:0000955",  # brain
            "UBERON:0002240",  # spinal cord
            "UBERON:0002048 (organoid)",  # lung (organoid)
        ]

        self.assertEqual(order_tissues(unordered_tissues), expected_ordered_tissues)
