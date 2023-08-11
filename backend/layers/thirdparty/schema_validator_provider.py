from typing import Tuple


class SchemaValidatorProviderInterface:
    def validate_and_save_labels(self, input_file: str, output_file: str) -> Tuple[bool, list, bool]:  # type: ignore
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
        from cellxgene_schema import validate

        return validate.validate(input_file, output_file)
