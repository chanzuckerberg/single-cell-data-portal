from dataclasses import asdict
from flask import g, make_response, jsonify

from backend.api_server.db import dbconnect
from backend.common.corpora_orm import DatasetArtifactFileType, DbCollection, IsPrimaryData, ProcessingStatus, ValidationStatus
from backend.common.utils.http_exceptions import (
    NotFoundHTTPException,
)
from backend.layers.api.router import get_business_logic
from backend.layers.api.transform import _dataset_asset_to_response, _dataset_processing_status_to_response
from backend.portal.api.app.v1.collections.collection_id.upload_links import upload_from_link
from backend.portal.api.collections_common import (
    get_dataset_else_error,
    delete_dataset_common,
)
from backend.portal.api.curation.v1.curation.collections.common import DATASET_ONTOLOGY_ELEMENTS, DATASET_ONTOLOGY_ELEMENTS_PREVIEW, EntityColumns, reshape_dataset_for_curation_api
from backend.layers.common.entities import CollectionId, DatasetArtifact, DatasetArtifactType, DatasetProcessingStatus, DatasetValidationStatus, DatasetVersion, DatasetVersionId


is_primary_data_mapping = {
    IsPrimaryData.PRIMARY: [True],
    IsPrimaryData.SECONDARY: [False],
    IsPrimaryData.BOTH: [True, False],
}


def _reshape_dataset_for_curation_api(d: DatasetVersion, preview=False) -> dict:
    # dataset = asdict(d)
    # if artifacts := dataset.pop("artifacts", []):
    #     dataset["dataset_assets"] = []
    #     for asset in artifacts:
    #         if asset["filetype"] in (DatasetArtifactFileType.H5AD, DatasetArtifactFileType.RDS):
    #             dataset["dataset_assets"].append(asset)
    # if processing_status := dataset.pop("processing_status", None):
    #     if processing_status["processing_status"] == ProcessingStatus.FAILURE:
    #         if processing_status["validation_status"] == ValidationStatus.INVALID:
    #             dataset["processing_status_detail"] = processing_status["validation_message"]
    #             dataset["processing_status"] = "VALIDATION_FAILURE"
    #         else:
    #             dataset["processing_status"] = "PIPELINE_FAILURE"
    #     else:
    #         dataset["processing_status"] = processing_status["processing_status"]
    # dataset_ontology_elements = DATASET_ONTOLOGY_ELEMENTS_PREVIEW if preview else DATASET_ONTOLOGY_ELEMENTS
    # for ontology_element in dataset_ontology_elements:
    #     if dataset_ontology_element := dataset.get(ontology_element):
    #         if not isinstance(dataset_ontology_element, list):
    #             # Package in array
    #             dataset[ontology_element] = [dataset_ontology_element]
    #     else:
    #         dataset[ontology_element] = []

    # if not preview:  # Add these fields only to full (and not preview) Dataset metadata response
    #     dataset["revision_of"] = dataset.pop("original_id", None)
    #     dataset["title"] = dataset.pop("name", None)
    #     if value := dataset.pop("is_primary_data", None):
    #         dataset["is_primary_data"] = is_primary_data_mapping.get(value, [])

    # return dataset

    artifacts = []
    for artifact in d.artifacts:
        if artifact.type in (DatasetArtifactType.H5AD, DatasetArtifactType.RDS):
            artifacts.append(_dataset_asset_to_response(artifact, d.dataset_id.id))



    dataset = {
        "assay": d.metadata.assay,
        "batch_condition": d.metadata.batch_condition,
        "cell_count": d.metadata.cell_count,
        "cell_type": d.metadata.cell_type,
        "dataset_assets": artifacts,
        "development_stage": d.metadata.development_stage,
        "disease": d.metadata.disease,
        "donor_id": d.metadata.donor_id,
        "explorer_url": "string",
        "id": d.dataset_id,
        "is_primary_data": d.metadata.is_primary_data,
        "mean_genes_per_cell": d.metadata.mean_genes_per_cell,
        "organism": d.metadata.organism,
        "processing_status": _dataset_processing_status_to_response(d.status, d.dataset_id.id),
        "processing_status_detail": "string",
        "revised_at": "string", # TODO
        "revision": 0,
        "schema_version": d.metadata.schema_version,
        "self_reported_ethnicity": d.metadata.self_reported_ethnicity,
        "sex": d.metadata.sex,
        "suspension_type": d.metadata.suspension_type,
        "tissue": d.metadata.tissue,
        "title": d.metadata.name,
        "x_approximate_distribution": d.metadata.x_approximate_distribution,
    }
            
    if d.status:
        if d.status.processing_status == DatasetProcessingStatus.FAILURE:
            if d.status.validation_status == DatasetValidationStatus.INVALID:
                dataset["processing_status_detail"] = d.status.validation_message
                dataset["processing_status"] = "VALIDATION_FAILURE"
            else:
                dataset["processing_status"] = "PIPELINE_FAILURE"

    dataset_ontology_elements = DATASET_ONTOLOGY_ELEMENTS_PREVIEW if preview else DATASET_ONTOLOGY_ELEMENTS
    for ontology_element in dataset_ontology_elements:
        if dataset_ontology_element := dataset.get(ontology_element):
            if not isinstance(dataset_ontology_element, list):
                # Package in array
                dataset[ontology_element] = [dataset_ontology_element]
        else:
            dataset[ontology_element] = []

    if not preview:  # Add these fields only to full (and not preview) Dataset metadata response
        dataset["revision_of"] = dataset.pop("original_id", None)
        dataset["title"] = dataset.pop("name", None)
        if value := dataset.pop("is_primary_data", None):
            dataset["is_primary_data"] = is_primary_data_mapping.get(value, [])

    return dataset


def get(collection_id: str, dataset_id: str = None):
    print("123")
    # business_logic = get_business_logic()

    # # TODO: double lookup
    # collection_version = business_logic.get_collection_version_from_canonical(CollectionId(collection_id))
    # if collection_version is None:
    #     raise NotFoundHTTPException("Collection not found!")

    # # TODO: double lookup
    # version = business_logic.get_dataset_version(DatasetVersionId(dataset_id))
    # if version is None:
    #     raise NotFoundHTTPException("Dataset not found")

    # response_body = _reshape_dataset_for_curation_api(version)
    # return make_response(jsonify(response_body), 200)


@dbconnect
def delete(token_info: dict, collection_id: str, dataset_id: str = None):
    db_session = g.db_session
    dataset = get_dataset_else_error(db_session, dataset_id, collection_id, include_tombstones=True)
    delete_dataset_common(db_session, dataset, token_info)
    return "", 202


def put(collection_id: str, dataset_id: str, body: dict, token_info: dict):
    upload_from_link(
        collection_id,
        token_info,
        body.get("url", body.get("link")),
        dataset_id,
    )
    return "", 202
