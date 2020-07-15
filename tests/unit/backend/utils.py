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
            name="Created Project",
            description="test",
            owner="test_user_id",
            s3_bucket="s3://fakebucket",
            tc_uri="https://fakeurl",
            needs_attestation=False,
            processing_state=ProcessingState.IN_VALIDATION.name,
            validation_state=ValidationState.NOT_VALIDATED.name,
        )
        bogus_data.update(**kwargs)
        return bogus_data
