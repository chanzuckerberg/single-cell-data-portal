from logging import Logger
from typing import List, Literal, Optional, Tuple
from backend.corpora.common.utils.corpora_constants import CorporaConstants
from backend.layers.business.business import BusinessLogic
from backend.layers.processing.prcess_logic import ProcessingLogic
from entities import DatasetConversionStatus, DatasetMetadata, DatasetProcessingStatus, DatasetStatusGeneric, DatasetStatusKey, DatasetUploadStatus, DatasetValidationStatus, DatasetVersionId, OntologyTermId

import scanpy
import numpy

class ProcessDownloadValidate(ProcessingLogic):

    def validate_h5ad_file_and_add_labels(self, dataset_id: DatasetVersionId, local_filename: str) -> Tuple[str, bool]:
        """
        Validates and labels the specified dataset file and updates the processing status in the database
        :param dataset_id: ID of the dataset to update
        :param local_filename: file name of the dataset to validate and label
        :return: file name of labeled dataset, boolean indicating if seurat conversion is possible
        """
        from cellxgene_schema import validate

        self.update_processing_status(dataset_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.VALIDATING)

        output_filename = CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME
        try:
            is_valid, errors, can_convert_to_seurat = validate.validate(local_filename, output_filename)
        except Exception as e:
            self.logger.error(f"Validation failed with exception: {e}!")
            status = dict(
                validation_status=ValidationStatus.INVALID,
                validation_message=str(e),
            )
            raise ValidationFailed(status) # TODO: why raise

        if not is_valid:
            self.logger.error(f"Validation failed with {len(errors)} errors!")
            status = dict(
                validation_status=ValidationStatus.INVALID,
                validation_message=errors,
            )
            raise ValidationFailed(status)
        else:
            self.logger.info("Validation complete")
            # TODO: optionally, these could be batched into one
            self.update_processing_status(dataset_id, DatasetStatusKey.H5AD, DatasetConversionStatus.CONVERTED)
            self.update_processing_status(dataset_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.VALID)
            return output_filename, can_convert_to_seurat


    def extract_metadata(self, filename) -> DatasetMetadata:
        """Pull metadata out of the AnnData file to insert into the dataset table."""

        adata = scanpy.read_h5ad(filename, backed="r")

        # TODO: Concern with respect to previous use of raising error when there is no raw layer.
        # This new way defaults to adata.X.
        if adata.raw is not None and adata.raw.X is not None:
            layer_for_mean_genes_per_cell = adata.raw.X
        else:
            layer_for_mean_genes_per_cell = adata.X

        # For mean_genes_per_cell, we only want the columns (genes) that have a feature_biotype of `gene`,
        # as opposed to `spike-in`
        filter_gene_vars = numpy.where(adata.var.feature_biotype == "gene")[0]

        # Calling np.count_nonzero on and h5py.Dataset appears to read the entire thing
        # into memory, so we need to chunk it to be safe.
        stride = 50000
        numerator, denominator = 0, 0
        for bounds in zip(
            range(0, layer_for_mean_genes_per_cell.shape[0], stride),
            range(stride, layer_for_mean_genes_per_cell.shape[0] + stride, stride),
        ):
            chunk = layer_for_mean_genes_per_cell[bounds[0] : bounds[1], filter_gene_vars]
            numerator += chunk.nnz if hasattr(chunk, "nnz") else numpy.count_nonzero(chunk)
            denominator += chunk.shape[0]

        def _get_term_pairs(base_term) -> List[OntologyTermId]:
            base_term_id = base_term + "_ontology_term_id"
            return [
                OntologyTermId(label=k[0], ontology_term_id=k[1])
                for k in adata.obs.groupby([base_term, base_term_id]).groups.keys()
            ]

        def _get_is_primary_data() -> Literal["PRIMARY", "SECONDARY", "BOTH"]:
            is_primary_data = adata.obs["is_primary_data"]
            if all(is_primary_data):
                return "PRIMARY"
            elif not any(is_primary_data):
                return "SECONDARY"
            else:
                return "BOTH"

        def _get_x_approximate_distribution() -> Optional[str]:
            if "X_approximate_distribution" in adata.uns:
                return adata.uns["X_approximate_distribution"].upper()
            else:
                return None

        def _get_batch_condition() -> Optional[str]:
            if "batch_condition" in adata.uns:
                return adata.uns["batch_condition"]
            else:
                return None

        return DatasetMetadata(
            name=adata.uns["title"],
            organism=_get_term_pairs("organism"),
            tissue=_get_term_pairs("tissue"),
            assay=_get_term_pairs("assay"),
            disease=_get_term_pairs("disease"),
            sex=_get_term_pairs("sex"),
            self_reported_ethnicity=_get_term_pairs("self_reported_ethnicity"),
            development_stage=_get_term_pairs("development_stage"),
            cell_count=adata.shape[0],
            mean_genes_per_cell=numerator / denominator,
            is_primary_data=_get_is_primary_data(),
            cell_type=_get_term_pairs("cell_type"),
            x_approximate_distribution=_get_x_approximate_distribution(), # TODO: pay attention
            schema_version=adata.uns["schema_version"],
            batch_condition=_get_batch_condition(), # TODO: pay attention
            donor_id=adata.obs["donor_id"].unique(),
            suspension_type=adata.obs["suspension_type"].unique(),
        )


    # TODO: not sure we need a method for this
    def wrapped_download_from_s3(self, dataset_id: DatasetVersionId, bucket_name: str, object_key: str, local_filename: str):
        """
        Wraps download_from_s3() to update the dataset's upload status
        :param dataset_id:
        :param bucket_name:
        :param object_key:
        :param local_filename:
        :return:
        """
        self.update_processing_status(dataset_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADING)
        self.download_from_s3(
            bucket_name=bucket_name,
            object_key=object_key,
            local_filename=local_filename,
        )
        self.update_processing_status(dataset_id, DatasetStatusKey.UPLOAD, DatasetUploadStatus.UPLOADED)


    def download_from_source_uri(self, dataset_id: DatasetVersionId, source_uri: str, local_path: str) -> str:
        """Given a source URI, download it to local_path.
        Handles fixing the url so it downloads directly.
        """

        file_url = self.uri_provider.parse(source_uri)
        if not file_url:
            raise ValueError(f"Malformed source URI: {source_uri}")

        # This is a bit ugly and should be done polymorphically instead, but Dropbox support will be dropped soon
        if file_url.scheme == "https":
            file_info = file_url.file_info()
            # TODO: unfortunately, downloader also needs a rewrite
            status = downloader.download(dataset_id, file_url.url, local_path, file_info["size"])
            self.logger.info(status)
        elif file_url.scheme == "s3":
            bucket_name = file_url.netloc
            key = self.remove_prefix(file_url.path, "/")
            self.wrapped_download_from_s3(
                dataset_id=dataset_id,
                bucket_name=bucket_name,
                object_key=key,
                local_filename=local_path,
            )
        else:
            raise ValueError(f"Download for URI scheme '{file_url.scheme}' not implemented")

        self.logger.info("Download complete") # TODO: remove
        return local_path


    # TODO: after upgrading to Python 3.9, replace this with removeprefix()
    def remove_prefix(self, string: str, prefix: str) -> str:
        if string.startswith(prefix):
            return string[len(prefix) :]
        else:
            return string[:]

    def process(self, dataset_id: DatasetVersionId, dropbox_url: str, artifact_bucket: str):
        """
        1. Download the original dataset from Dropbox
        2. Validate and label it
        3. Upload the labeled dataset to the artifact bucket
        :param dataset_id:
        :param dropbox_url:
        :param artifact_bucket:
        :return:
        """

        self.update_processing_status(dataset_id, DatasetStatusKey.PROCESSING, DatasetProcessingStatus.PENDING)

        # Download the original dataset from Dropbox
        local_filename = self.download_from_source_uri(
            dataset_id=dataset_id,
            source_uri=dropbox_url,
            local_path=CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME,
        )

        # Validate and label the dataset
        file_with_labels, can_convert_to_seurat = self.validate_h5ad_file_and_add_labels(dataset_id, local_filename)
        # Process metadata
        metadata = self.extract_metadata(file_with_labels)
        self.business_logic.set_dataset_metadata(dataset_id, metadata)

        if not can_convert_to_seurat:
            self.update_processing_status(dataset_id, DatasetStatusKey.RDS, DatasetConversionStatus.SKIPPED)
            self.logger.info(f"Skipping Seurat conversion for dataset {dataset_id}")

        bucket_prefix = self.get_bucket_prefix(dataset_id)
        # Upload the original dataset to the artifact bucket
        self.create_artifact(
            local_filename, "H5AD", bucket_prefix, dataset_id, artifact_bucket, DatasetStatusKey.H5AD
        )
        # Upload the labeled dataset to the artifact bucket
        self.create_artifact(
            file_with_labels, "H5AD", bucket_prefix, dataset_id, artifact_bucket, DatasetStatusKey.H5AD
        )
