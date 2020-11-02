from backend.corpora.common.corpora_orm import (
    CollectionVisibility,
)


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
            organism={"ontology_term_id": "123", "label": "organism"},
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
            sex=["male", "female", "mixed"],
            ethnicity=[{"ontology_term_id": "", "label": "unknown"}],
            development_stage=[{"ontology_term_id": "HsapDv:0011", "label": "just a baby"}],
            collection_id="test_collection_id",
            collection_visibility=CollectionVisibility.PUBLIC.name,
        )

        bogus_data.update(**kwargs)
        return bogus_data
