from typing import List, Protocol, Tuple

from cellxgene_schema import validate
from cellxgene_schema.migrate import migrate
from cellxgene_schema.schema import get_current_schema_version


class SchemaValidatorProviderInterface(Protocol):
    def validate_and_save_labels(self, input_file: str, output_file: str) -> Tuple[bool, list, bool]:
        pass

    def migrate(self, input_file: str, output_file: str, collection_id: str, dataset_id: str) -> List[str]:
        """
        Runs `cellxgene-schema migrate` on the provided `input_file`.
        """
        pass

    def get_current_schema_version(self) -> str:
        """
        Returns the current schema version
        """
        pass


class SchemaValidatorProvider(SchemaValidatorProviderInterface):
    def validate_and_save_labels(self, input_file: str, output_file: str) -> Tuple[bool, list, bool]:
        """
        Runs `cellxgene-schema validate` on the provided `input_file`. This also saves a labeled copy
        of the artifact to `output_file`.
        Returns a tuple that contains, in order:
        1. A boolean that indicates whether the artifact is valid
        2. A List[str] with the validation errors. This is only defined if the first boolean is false
        3. A boolean that indicates whether the artifact is Seurat convertible
        """

        return validate.validate(input_file, output_file)

    def migrate(self, input_file, output_file, collection_id, dataset_id) -> List[str]:
        """
        Runs `cellxgene-schema migrate` on the provided `input_file`.
        """
        return migrate(input_file, output_file, collection_id, dataset_id)

    def get_current_schema_version(self) -> str:
        return get_current_schema_version()
