from typing import List, Literal, Optional, Tuple

import numpy
import scanpy

from backend.common.corpora_config import CorporaConfig
from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import (
    CollectionVersionId,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetMetadata,
    DatasetProcessingStatus,
    DatasetStatusKey,
    DatasetUploadStatus,
    DatasetValidationStatus,
    DatasetVersionId,
    OntologyTermId,
)
from backend.layers.processing.downloader import Downloader
from backend.layers.processing.exceptions import ValidationFailed
from backend.layers.processing.logger import logit
from backend.layers.processing.process_logic import ProcessingLogic
from backend.layers.thirdparty.s3_provider import S3ProviderInterface
from backend.layers.thirdparty.schema_validator_provider import SchemaValidatorProviderInterface
from backend.layers.thirdparty.uri_provider import UriProviderInterface


class ProcessDownloadValidate(ProcessingLogic):
    """
    Base class for handling the `Download and Validate` step of the step function.
    This will:
    1. Download the original artifact from the provided URI
    2. Run the cellxgene-schema validator
    3. Save and upload a labeled copy of the original artifact (local.h5ad)
    4. Upload a copy of the original artifact (raw.h5ad)
    5. Persist the dataset metadata on the database
    6. Determine if a Seurat conversion is possible (it is not if the X matrix has more than 2**32-1 nonzero values)
    If this step completes successfully, ProcessCxg and ProcessSeurat will start in parallel.
    If this step fails, the handle_failures lambda will be invoked.
    """

    downloader: Downloader
    schema_validator: SchemaValidatorProviderInterface

    def __init__(
        self,
        business_logic: BusinessLogicInterface,
        uri_provider: UriProviderInterface,
        s3_provider: S3ProviderInterface,
        downloader: Downloader,
        schema_validator: SchemaValidatorProviderInterface,
    ) -> None:
        super().__init__()
        self.business_logic = business_logic
        self.uri_provider = uri_provider
        self.s3_provider = s3_provider
        self.downloader = downloader
        self.schema_validator = schema_validator

    @logit
    def validate_h5ad_file_and_add_labels(
        self, collection_id: CollectionVersionId, dataset_id: DatasetVersionId, local_filename: str
    ) -> Tuple[str, bool]:
        """
        Validates and labels the specified dataset file and updates the processing status in the database
        :param dataset_id: ID of the dataset to update
        :param local_filename: file name of the dataset to validate and label
        :return: file name of labeled dataset, boolean indicating if seurat conversion is possible
        """
        # TODO: use a provider here

        self.update_processing_status(dataset_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.VALIDATING)

        output_filename = CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME
        try:
            is_valid, errors, can_convert_to_seurat = self.schema_validator.validate_and_save_labels(
                local_filename, output_filename
            )
        except Exception as e:
            self.logger.exception("validation failed")
            raise ValidationFailed([str(e)]) from None

        if not is_valid:
            raise ValidationFailed(errors)
        else:
            if CorporaConfig().schema_4_feature_flag.lower() == "true":
                self.populate_dataset_citation(collection_id, dataset_id, output_filename)

            # TODO: optionally, these could be batched into one
            self.update_processing_status(dataset_id, DatasetStatusKey.H5AD, DatasetConversionStatus.CONVERTED)
            self.update_processing_status(dataset_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.VALID)
            return output_filename, can_convert_to_seurat

    def populate_dataset_citation(
        self, collection_id: CollectionVersionId, dataset_id: DatasetVersionId, adata_path: str
    ) -> None:
        """
        Builds citation string and updates the 'uns' dict of the adata at adata_path

        :param collection_id: version ID for collection dataset is being uploaded to
        :param dataset_id: version ID for dataset
        :param adata_path: filepath to adata object that will be updated with citation
        """
        dataset_assets_base_url = CorporaConfig().dataset_assets_base_url
        collections_base_url = CorporaConfig().collections_base_url
        citation = ""
        collection = self.business_logic.get_collection_version(collection_id)
        doi = next((link.uri for link in collection.metadata.links if link.type == "DOI"), None)
        if doi:
            citation += f"Publication: {doi} "
        citation += f"Dataset Version: {dataset_assets_base_url}/{dataset_id}.h5ad "
        citation += (
            f"curated and distributed by CZ CELLxGENE Discover in Collection: "
            f"{collections_base_url}/{collection_id}"
        )
        adata = scanpy.read_h5ad(adata_path)
        adata.uns["citation"] = citation
        adata.write(adata_path)

    @logit
    def extract_metadata(self, filename) -> DatasetMetadata:
        """Pull metadata out of the AnnData file to insert into the dataset table."""

        adata = scanpy.read_h5ad(filename, backed="r")

        # TODO: Concern with respect to previous use of raising error when there is no raw layer.
        # This new way defaults to adata.X.
        layer_for_mean_genes_per_cell = adata.raw.X if adata.raw is not None and adata.raw.X is not None else adata.X

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
                for k in adata.obs.groupby([base_term, base_term_id]).groups
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
            primary_cell_count=int(adata.obs["is_primary_data"].astype("int").sum()),
            mean_genes_per_cell=numerator / denominator,
            is_primary_data=_get_is_primary_data(),
            cell_type=_get_term_pairs("cell_type"),
            x_approximate_distribution=_get_x_approximate_distribution(),
            schema_version=adata.uns["schema_version"],
            batch_condition=_get_batch_condition(),
            donor_id=adata.obs["donor_id"].unique(),
            suspension_type=adata.obs["suspension_type"].unique(),
            feature_count=adata.var.shape[0],
            feature_biotype=adata.var["feature_biotype"].unique(),
            feature_reference=adata.var["feature_reference"].unique(),
            default_embedding=adata.uns.get("default_embedding"),
            embeddings=adata.obsm_keys(),
            raw_data_location="raw.X" if adata.raw else "X",
            citation=adata.uns.get("citation"),
        )

    def wrapped_download_from_s3(
        self, dataset_id: DatasetVersionId, bucket_name: str, object_key: str, local_filename: str
    ):
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

    @logit
    def download_from_source_uri(self, dataset_id: DatasetVersionId, source_uri: str, local_path: str) -> str:
        """Given a source URI, download it to local_path.
        Handles fixing the url so it downloads directly.
        """
        self.logger.info("Start download")
        file_url = self.uri_provider.parse(source_uri)
        if not file_url:
            raise ValueError(f"Malformed source URI: {source_uri}")

        # This is a bit ugly and should be done polymorphically instead, but Dropbox support will be dropped soon
        if file_url.scheme == "https":
            file_info = self.uri_provider.get_file_info(source_uri)
            status = self.downloader.download(dataset_id, file_url.url, local_path, file_info.size)
            self.logger.info(status)  # TODO: this log is awful
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

        self.logger.info("Download complete")  # TODO: remove
        return local_path

    # TODO: after upgrading to Python 3.9, replace this with removeprefix()
    def remove_prefix(self, string: str, prefix: str) -> str:
        if string.startswith(prefix):
            return string[len(prefix) :]
        else:
            return string[:]

    def process(
        self,
        collection_id: CollectionVersionId,
        dataset_id: DatasetVersionId,
        dropbox_url: str,
        artifact_bucket: str,
        datasets_bucket: str,
    ):
        """
        1. Download the original dataset from Dropbox
        2. Validate and label it
        3. Upload the labeled dataset to the artifact bucket
        4. Upload the labeled dataset to the datasets bucket
        :param collection_id
        :param dataset_id:
        :param dropbox_url:
        :param artifact_bucket:
        :param datasets_bucket:
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
        file_with_labels, can_convert_to_seurat = self.validate_h5ad_file_and_add_labels(
            collection_id, dataset_id, local_filename
        )
        # Process metadata
        metadata = self.extract_metadata(file_with_labels)
        self.business_logic.set_dataset_metadata(dataset_id, metadata)

        if not can_convert_to_seurat:
            self.update_processing_status(dataset_id, DatasetStatusKey.RDS, DatasetConversionStatus.SKIPPED)
            self.logger.info(f"Skipping Seurat conversion for dataset {dataset_id}")

        key_prefix = self.get_key_prefix(dataset_id.id)
        # Upload the original dataset to the artifact bucket
        self.create_artifact(
            local_filename,
            DatasetArtifactType.RAW_H5AD,
            key_prefix,
            dataset_id,
            artifact_bucket,
            DatasetStatusKey.H5AD,
        )
        # Upload the labeled dataset to the artifact bucket
        self.create_artifact(
            file_with_labels,
            DatasetArtifactType.H5AD,
            key_prefix,
            dataset_id,
            artifact_bucket,
            DatasetStatusKey.H5AD,
            datasets_bucket=datasets_bucket,
        )
