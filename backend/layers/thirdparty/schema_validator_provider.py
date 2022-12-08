from typing import Tuple


class SchemaValidatorProviderInterface:
    def validate_and_save_labels(self, input_file: str, output_file: str) -> Tuple[bool, list, bool]:
        pass


class SchemaValidatorProvider(SchemaValidatorProviderInterface):
    def validate_and_save_labels(self, input_file: str, output_file: str) -> Tuple[bool, list, bool]:
        from cellxgene_schema import validate

        return validate.validate(input_file, output_file)
