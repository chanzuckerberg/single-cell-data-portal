from backend.corpora.common.corpora_orm import (
    ProcessingState,
    ValidationState,
    ProjectStatus,
)


class BogusProjectParams:
    @classmethod
    def get(cls, **kwargs):
        bogus_data = dict(
            status=ProjectStatus.EDIT.name,
            owner="test_user_id",
            processing_state=ProcessingState.IN_VALIDATION.name,
            validation_state=ValidationState.NOT_VALIDATED.name,
        )
        bogus_data.update(**kwargs)
        return bogus_data


class BogusDatasetParams:
    @classmethod
    def get(cls, **kwargs):
        bogus_data = dict(
            name="create_dataset",
            organism={"ontology_term_id": "123", "label": "organism"},
            tissue=[{"ontology_term_id": "UBERON:1111", "label": "brain"},
                    {"ontology_term_id": "UBERON:2222", "label": "something unusual"}],
            assay=[{"ontology_term_id": "ABC:00000123", "label": "10x"}],
            disease=[{"ontology_term_id": "MONDO:0000456"}, {"label": "heart disease"},
                     {"ontology_term_id": "MONDO:0000789"}, {"label": "lung disease"}],
            sex=["male", "female", "mixed"],
            ethnicity=[{"ontology_term_id": "", "label": "unknown"}],
            development_stage=[{"ontology_term_id": "HsapDv:0011", "label": "just a baby"}],
            source_data_location="location",
            preprint_doi="preprint",
            publication_doi="publication",
            project_id="test_project_id",
            project_status=ProjectStatus.LIVE.name,
        )

        bogus_data.update(**kwargs)
        return bogus_data
