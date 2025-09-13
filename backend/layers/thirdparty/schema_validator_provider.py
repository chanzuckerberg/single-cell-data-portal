from typing import List, Optional, Protocol, Tuple

from cellxgene_schema import validate
from cellxgene_schema.migrate import migrate
from cellxgene_schema.schema import get_current_schema_version

from backend.layers.processing.exceptions import AddLabelsFailed


class SchemaValidatorProviderInterface(Protocol):
    def validate_anndata(self, input_file: str) -> Tuple[bool, list, bool]:
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

    def add_labels(self, input_file: str, output_file: str) -> None:
        """
        Adds labels to the provided `input_file` and writes the result to `output_file`.
        """
        pass

    def validate_atac(self, fragment_file, anndata_file, output_file) -> Tuple[Optional[List[str]], str, str]:
        """
        Validates an ATAC fragment file against an anndata file.

        Returns a tuple that contains, in order:
        1. A List[str] with the validation errors. This is only defined if the first boolean is false
        2. The path to the index file
        3. The path to the fragment file
        """
        pass

    def check_anndata_requires_fragment(self, anndata_file) -> bool:
        """
        Check if an anndata file requires a fragment file
        """
        pass

    def deduplicate_fragments(self, fragment_file: str) -> str:
        """
        Deduplicates and replaces the provided `fragment_file`.
        """
        pass


class SchemaValidatorProvider(SchemaValidatorProviderInterface):
    def validate_anndata(self, input_file: str) -> Tuple[bool, list, bool]:
        """
        Runs `cellxgene-schema validate` on the provided `input_file`.
        Returns a tuple that contains, in order:
        1. A boolean that indicates whether the artifact is valid
        2. A List[str] with the validation errors. This is only defined if the first boolean is false
        3. A boolean that indicates whether the artifact is Seurat convertible
        """

        return validate.validate(input_file)

    def migrate(self, input_file, output_file, collection_id, dataset_id) -> List[str]:
        """
        Runs `cellxgene-schema migrate` on the provided `input_file`.
        """
        return migrate(input_file, output_file, collection_id, dataset_id)

    def get_current_schema_version(self) -> str:
        return get_current_schema_version()

    def add_labels(self, input_file: str, output_file: str) -> None:
        """
        Adds labels to the provided `input_file` and writes the result to `output_file`.
        """
        from cellxgene_schema.utils import read_h5ad
        from cellxgene_schema.write_labels import AnnDataLabelAppender

        adata = read_h5ad(input_file)
        anndata_label_adder = AnnDataLabelAppender(adata)
        if not anndata_label_adder.write_labels(output_file):
            raise AddLabelsFailed(anndata_label_adder.errors)

    def validate_atac(self, fragment_file, anndata_file, output_file) -> Tuple[Optional[List[str]], str, str]:
        """
        Validates an ATAC fragment file against an anndata file.

        Returns a tuple that contains, in order:
        1. A List[str] with the validation errors. This is only defined if the first boolean is false
        2. The path to the index file
        3. The path to the fragment file
        """
        import cellxgene_schema.atac_seq as atac_seq

        index_file = output_file + ".tbi"
        return (
            atac_seq.process_fragment(fragment_file, anndata_file, True, output_file=output_file),
            index_file,
            output_file,
        )

    def check_anndata_requires_fragment(self, anndata_file) -> bool:
        import cellxgene_schema.atac_seq as atac_seq

        return atac_seq.check_anndata_requires_fragment(anndata_file)

    def deduplicate_fragments(self, fragment_file: str) -> str:
        import os

        import cellxgene_schema.atac_seq as atac_seq

        output_file = atac_seq.deduplicate_fragment_rows(fragment_file)
        if not output_file or not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
            raise RuntimeError("Deduplication failed: output file not created. Original file not removed.")
        os.remove(fragment_file)
        return output_file
