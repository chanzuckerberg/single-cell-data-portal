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
            organism="organism",
            organism_ontology="123",
            tissue="tissue",
            tissue_ontology="123",
            assay="assay",
            assay_ontology="123",
            disease="diseas",
            disease_ontology="123",
            sex="F",
            ethnicity="ethnicity",
            ethnicity_ontology="123",
            source_data_location="location",
            preprint_doi="preprint",
            publication_doi="publication",
            project_id="test_project_id",
            project_status=ProjectStatus.LIVE.name,
        )

        bogus_data.update(**kwargs)
        return bogus_data
