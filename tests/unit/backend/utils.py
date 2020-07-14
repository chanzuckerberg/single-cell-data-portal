from backend.corpora.common.corpora_orm import (
    ProcessingState,
    ValidationState,
    ProjectStatus,
)


class ProjectParams:
    @classmethod
    def get(cls):
        return dict(
            status=ProjectStatus.EDIT.name,
            name="Created Project",
            description="test",
            owner="test_user_id",
            s3_bucket="s3://fakebucket",
            tc_uri="https://fakeurl",
            needs_attestation=False,
            processing_state=ProcessingState.IN_VALIDATION.name,
            validation_state=ValidationState.NOT_VALIDATED.name,
        )
