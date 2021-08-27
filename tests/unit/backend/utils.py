import random
import string

from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
    ConversionStatus,
    IsPrimaryData,
    UploadStatus,
    ValidationStatus,
    XApproximateDistribution,
)


class BogusProcessingStatusParams:
    @classmethod
    def get(cls, **kwargs):
        bogus_data = dict(
            upload_status=UploadStatus.UPLOADING,
            upload_progress=1 / 9,
            validation_status=ValidationStatus.NA,
            conversion_loom_status=ConversionStatus.NA,
            conversion_rds_status=ConversionStatus.NA,
            conversion_cxg_status=ConversionStatus.NA,
            conversion_anndata_status=ConversionStatus.NA,
        )
        bogus_data.update(**kwargs)
        return bogus_data


class BogusCollectionParams:
    @classmethod
    def get(cls, **kwargs):
        bogus_data = dict(
            visibility=CollectionVisibility.PRIVATE.name, owner="test_user_id", data_submission_policy_version=0
        )
        bogus_data.update(**kwargs)
        return bogus_data


class BogusDatasetParams:
    @classmethod
    def get(cls, **kwargs):
        bogus_data = dict(
            name="create_dataset",
            organism=[
                {"ontology_term_id": "123", "label": "organism"}
            ],
            tissue=[
                {"ontology_term_id": "UBERON:1111", "label": "brain"},
                {"ontology_term_id": "UBERON:2222", "label": "something unusual"},
            ],
            assay=[{"ontology_term_id": "ABC:00000123", "label": "10x"}],
            disease=[
                {"ontology_term_id": "MONDO:0000456"},
                {"label": "heart disease"},
                {"ontology_term_id": "MONDO:0000789"},
                {"label": "lung disease"},
            ],
            sex=[
                {"ontology_term_id": "M", "label": "male"},
                {"ontology_term_id": "F", "label": "female"},
                {"ontology_term_id": "MF", "label": "mixed"},
            ],
            ethnicity=[{"ontology_term_id": "", "label": "unknown"}],
            development_stage=[{"ontology_term_id": "HsapDv:0011", "label": "just a baby"}],
            cell_type=[{"ontology_term_id": "Hepatic-1A", "label": "liver"}],
            is_primary_data=IsPrimaryData.PRIMARY.name,
            collection_id="test_collection_id",
            collection_visibility=CollectionVisibility.PUBLIC.name,
            explorer_url="test_url",
            x_normalization="normal",
            x_approximate_distribution=XApproximateDistribution.NORMAL.name,
            schema_version="2.0.0",
        )

        bogus_data.update(**kwargs)
        return bogus_data


class BogusGenesetParams:
    @classmethod
    def get(cls, gene_count=6, **kwargs):
        genes = []
        for i in range(gene_count):
            gene = {
                "gene_symbol": f"{i}",
                "gene_description": "describe a gene",
                "additional_params": {},
            }
            if i % 3 == 0:
                gene["additional_params"] = {
                    "provenance1": "some words",
                    "provenance1_description": "another set of words",
                }
            genes.append(gene)
        bogus_data = dict(
            description="This is a geneset bwhahaha",
            name=cls.generate_random_string(7),
            genes=genes,
            collection_id="test_collection_id",
            collection_visibility=CollectionVisibility.PUBLIC.name,
        )

        bogus_data.update(**kwargs)
        return bogus_data

    @staticmethod
    def generate_random_string(length=7):
        return "".join(random.choice(string.ascii_letters) for i in range(length))
