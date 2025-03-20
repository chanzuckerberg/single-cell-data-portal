from typing import Any, Dict, List, Literal, Optional

import h5py
import numpy
from cellxgene_schema.utils import read_h5ad

from backend.common.utils.corpora_constants import CorporaConstants
from backend.layers.business.business_interface import BusinessLogicInterface
from backend.layers.common.entities import (
    CollectionVersionId,
    DatasetArtifactType,
    DatasetConversionStatus,
    DatasetMetadata,
    DatasetStatusKey,
    DatasetValidationStatus,
    DatasetVersionId,
    OntologyTermId,
    SpatialMetadata,
    TissueOntologyTermId,
)
from backend.layers.processing.exceptions import AddLabelsFailed
from backend.layers.processing.logger import logit
from backend.layers.processing.process_logic import ProcessingLogic
from backend.layers.thirdparty.s3_provider import S3ProviderInterface
from backend.layers.thirdparty.schema_validator_provider import SchemaValidatorProviderInterface
from backend.layers.thirdparty.uri_provider import UriProviderInterface


class ProcessAddLabels(ProcessingLogic):
    """
    Base class for handling the `add label` step of the step function.
    This will:
        1. Download the h5ad artifact
        2. Add labels to h5ad using cellxgene-schema
        3. Persist the dataset metadata on the database
        4. upload the labeled file to S3
        5. set DatasetStatusKey.H5AD status to DatasetUploadStatus.UPLOADED


    If this step completes successfully, ProcessCxg and ProcessSeurat will start in parallel.
    If this step fails, the handle_failures lambda will be invoked.
    """

    schema_validator: SchemaValidatorProviderInterface

    def __init__(
        self,
        business_logic: BusinessLogicInterface,
        uri_provider: UriProviderInterface,
        s3_provider: S3ProviderInterface,
        schema_validator: SchemaValidatorProviderInterface,
    ) -> None:
        super().__init__()
        self.business_logic = business_logic
        self.uri_provider = uri_provider
        self.s3_provider = s3_provider
        self.schema_validator = schema_validator

    @logit
    def add_labels(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId, local_filename: str
    ) -> str:
        """
        labels the specified dataset file and updates the processing status in the database
        :param dataset_version_id: version ID of the dataset to update
        :param collection_version_id: version ID of the collection dataset is being uploaded to
        :param local_filename: file name of the dataset to validate and label
        :return: file name of labeled dataset, boolean indicating if seurat conversion is possible
        """
        output_filename = self.get_file_path(CorporaConstants.LABELED_H5AD_ARTIFACT_FILENAME)
        self.schema_validator.add_labels(local_filename, output_filename)
        self.populate_dataset_citation(collection_version_id, dataset_version_id, output_filename)
        self.update_processing_status(dataset_version_id, DatasetStatusKey.H5AD, DatasetConversionStatus.CONVERTED)
        return output_filename

    def populate_dataset_citation(
        self, collection_version_id: CollectionVersionId, dataset_version_id: DatasetVersionId, adata_path: str
    ) -> None:
        """
        Builds citation string and updates the 'uns' dict of the adata at adata_path

        :param collection_version_id: version ID for collection dataset is being uploaded to
        :param dataset_version_id: version ID for dataset
        :param adata_path: filepath to adata object that will be updated with citation
        """
        collection = self.business_logic.get_collection_version(collection_version_id, get_tombstoned=False)
        doi = next((link.uri for link in collection.metadata.links if link.type == "DOI"), None)
        citation = self.business_logic.generate_dataset_citation(collection.collection_id, dataset_version_id, doi)
        with h5py.File(adata_path, "r+") as f:
            f["uns"].create_dataset("citation", data=citation)

    def get_spatial_metadata(self, spatial_dict: Dict[str, Any]) -> Optional[SpatialMetadata]:
        """
        Extracts spatial dataset metadata from the uns dict of an AnnData object

        :param spatial_dict: the value of the 'spatial' key from the uns dict of an AnnData object
        :return: SpatialMetadata object
        """
        is_single = spatial_dict.get("is_single")
        has_fullres = False
        spatial_library_ids = [key for key in spatial_dict if key != "is_single"]
        # schema validation ensures there can only be at max, one other key in uns["spatial"] if "is_single" is True
        library_id = spatial_library_ids.pop() if spatial_library_ids else None
        if library_id and "images" in spatial_dict[library_id] and "fullres" in spatial_dict[library_id]["images"]:
            has_fullres = True
        return SpatialMetadata(is_single=bool(is_single), has_fullres=has_fullres)

    @logit
    def extract_metadata(self, filename) -> DatasetMetadata:
        """Pull metadata out of the AnnData file to insert into the dataset table."""

        adata = read_h5ad(filename)

        layer_for_mean_genes_per_cell = adata.raw.X if adata.raw is not None and adata.raw.X is not None else adata.X

        # For mean_genes_per_cell, we only want the columns (genes) that have a feature_biotype of `gene`,
        # as opposed to `spike-in`
        filter_gene_vars = numpy.where(adata.var.feature_biotype == "gene")[0]
        filtered_matrix = layer_for_mean_genes_per_cell[:, filter_gene_vars]
        nnz_gene_exp = self.schema_validator.count_matrix_nonzero(filtered_matrix)
        total_cells = layer_for_mean_genes_per_cell.shape[0]
        mean_genes_per_cell = nnz_gene_exp / total_cells

        def _get_term_pairs(base_term) -> List[OntologyTermId]:
            base_term_id = base_term + "_ontology_term_id"
            return [
                OntologyTermId(label=k[0], ontology_term_id=k[1])
                for k in adata.obs.groupby([base_term, base_term_id]).groups
            ]

        def _get_tissue_terms() -> List[TissueOntologyTermId]:
            return [
                TissueOntologyTermId(label=k[0], ontology_term_id=k[1], tissue_type=k[2])
                for k in adata.obs.groupby(["tissue", "tissue_ontology_term_id", "tissue_type"]).groups
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
            tissue=_get_tissue_terms(),
            assay=_get_term_pairs("assay"),
            disease=_get_term_pairs("disease"),
            sex=_get_term_pairs("sex"),
            self_reported_ethnicity=_get_term_pairs("self_reported_ethnicity"),
            development_stage=_get_term_pairs("development_stage"),
            cell_count=adata.shape[0],
            primary_cell_count=int(adata.obs["is_primary_data"].astype("int").sum()),
            mean_genes_per_cell=mean_genes_per_cell,
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
            spatial=self.get_spatial_metadata(adata.uns["spatial"]) if "spatial" in adata.uns else None,
        )

    def process(
        self,
        collection_version_id: CollectionVersionId,
        dataset_version_id: DatasetVersionId,
        artifact_bucket: str,
        datasets_bucket: str,
    ):
        """
        1. Download the original dataset from URI
        2. Validate and label it
        3. Upload the labeled dataset to the artifact bucket
        4. Upload the labeled dataset to the datasets bucket
        :param collection_version_id
        :param dataset_version_id:
        :param artifact_bucket:
        :param datasets_bucket:
        :return:
        """
        self.update_processing_status(dataset_version_id, DatasetStatusKey.VALIDATION, DatasetValidationStatus.VALID)
        # Download the original dataset from S3
        key_prefix = self.get_key_prefix(dataset_version_id.id)
        original_h5ad_artifact_file_name = self.get_file_path(CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME)
        object_key = f"{key_prefix}/{CorporaConstants.ORIGINAL_H5AD_ARTIFACT_FILENAME}"
        self.download_from_s3(artifact_bucket, object_key, original_h5ad_artifact_file_name)

        # label the dataset
        try:
            file_with_labels = self.add_labels(
                collection_version_id, dataset_version_id, original_h5ad_artifact_file_name
            )
        except Exception as e:
            self.logger.exception(f"An unexpected error occurred while adding labels to the data set: {e}")
            raise AddLabelsFailed() from e
        # Process metadata
        metadata = self.extract_metadata(file_with_labels)
        self.business_logic.set_dataset_metadata(dataset_version_id, metadata)
        # Upload the labeled dataset to the artifact bucket
        self.create_artifact(
            file_with_labels,
            DatasetArtifactType.H5AD,
            key_prefix,
            dataset_version_id,
            DatasetStatusKey.H5AD,
            artifact_bucket,
            datasets_bucket=datasets_bucket,
        )
