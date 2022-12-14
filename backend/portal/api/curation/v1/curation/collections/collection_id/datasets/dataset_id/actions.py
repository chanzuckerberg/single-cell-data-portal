from flask import Response, jsonify, make_response
from backend.common.corpora_orm import IsPrimaryData
from backend.common.utils.exceptions import MaxFileSizeExceededException
from backend.common.utils.http_exceptions import ForbiddenHTTPException, InvalidParametersHTTPException, MethodNotAllowedException, NotFoundHTTPException, TooLargeHTTPException
from backend.layers.api.router import get_business_logic
from backend.layers.api.transform import dataset_asset_to_response, dataset_processing_status_to_response, ontology_term_ids_to_response
from backend.layers.auth.user_info import UserInfo
from backend.layers.business.exceptions import CollectionIsPublishedException, CollectionNotFoundException, CollectionUpdateException, DatasetInWrongStatusException, DatasetNotFoundException, InvalidURIException

from backend.layers.common.entities import CollectionId, CollectionVersionId, DatasetArtifact, DatasetArtifactType, DatasetProcessingStatus, DatasetValidationStatus, DatasetVersion, DatasetVersionId
from backend.portal.api.curation.v1.curation.collections.common import DATASET_ONTOLOGY_ELEMENTS, DATASET_ONTOLOGY_ELEMENTS_PREVIEW, get_infered_dataset_version


is_primary_data_mapping = {
    "PRIMARY": [True],
    "SECONDARY": [False],
    "BOTH": [True, False],
}


def _reshape_dataset_for_curation_api(d: DatasetVersion, preview=False) -> dict:
    artifacts = []
    for artifact in d.artifacts:
        if artifact.type in (DatasetArtifactType.H5AD, DatasetArtifactType.RDS):
            artifacts.append(dataset_asset_to_response(artifact, d.dataset_id.id))

    dataset = {
        "assay": ontology_term_ids_to_response(d.metadata.assay),
        "batch_condition": d.metadata.batch_condition,
        "cell_count": d.metadata.cell_count,
        "cell_type": ontology_term_ids_to_response(d.metadata.cell_type),
        "dataset_assets": artifacts,
        "development_stage": ontology_term_ids_to_response(d.metadata.development_stage),
        "disease": ontology_term_ids_to_response(d.metadata.disease),
        "donor_id": d.metadata.donor_id,
        "explorer_url": "string",
        "id": d.dataset_id.id,
        "mean_genes_per_cell": d.metadata.mean_genes_per_cell,
        "organism": ontology_term_ids_to_response(d.metadata.organism),
        "processing_status": dataset_processing_status_to_response(d.status, d.dataset_id.id),
        "processing_status_detail": "string",
        "revised_at": "string", # TODO
        "revision": 0,
        "schema_version": d.metadata.schema_version,
        "self_reported_ethnicity": ontology_term_ids_to_response(d.metadata.self_reported_ethnicity),
        "sex": ontology_term_ids_to_response(d.metadata.sex),
        "suspension_type": d.metadata.suspension_type,
        "tissue": ontology_term_ids_to_response(d.metadata.tissue),
        "x_approximate_distribution": d.metadata.x_approximate_distribution.upper(),
    }

    if d.status:
        if d.status.processing_status == DatasetProcessingStatus.FAILURE:
            if d.status.validation_status == DatasetValidationStatus.INVALID:
                dataset["processing_status_detail"] = d.status.validation_message
                dataset["processing_status"] = "VALIDATION_FAILURE"
            else:
                dataset["processing_status"] = "PIPELINE_FAILURE"
        else:
            dataset["processing_status"] = d.status.processing_status

    dataset_ontology_elements = DATASET_ONTOLOGY_ELEMENTS_PREVIEW if preview else DATASET_ONTOLOGY_ELEMENTS
    for ontology_element in dataset_ontology_elements:
        if dataset_ontology_element := dataset.get(ontology_element):
            if not isinstance(dataset_ontology_element, list):
                # Package in array
                dataset[ontology_element] = [dataset_ontology_element]
        else:
            dataset[ontology_element] = []

    if not preview:  # Add these fields only to full (and not preview) Dataset metadata response
        dataset["revision_of"] = d.canonical_dataset.dataset_id.id
        dataset["title"] = d.metadata.name
        if d.metadata.is_primary_data is not None:
            dataset["is_primary_data"] = is_primary_data_mapping.get(d.metadata.is_primary_data, [])

    return dataset


def get(collection_id: str, dataset_id: str = None):
    business_logic = get_business_logic()

    collection_version = business_logic.get_collection_version_from_canonical(CollectionId(collection_id))
    if collection_version is None:
        raise NotFoundHTTPException("Collection not found!")

    version = get_infered_dataset_version(dataset_id)
    if version is None:
        raise NotFoundHTTPException("Dataset not found")

    response_body = _reshape_dataset_for_curation_api(version)
    return make_response(jsonify(response_body), 200)


def delete(token_info: dict, collection_id: str, dataset_id: str = None):

    business_logic = get_business_logic()
    user_info = UserInfo(token_info)

    # TODO: deduplicate from ApiCommon. We need to settle the class/module level debate before can do that
    dataset_version = business_logic.get_dataset_version(DatasetVersionId(dataset_id))
    if dataset_version is None:
        raise ForbiddenHTTPException()

    collection_version = business_logic.get_collection_version(CollectionVersionId(collection_id))
    # If the collection does not exist, it means that the dataset is orphaned and therefore we cannot
    # determine the owner. This should not be a problem - we won't need its state at that stage.
    if collection_version is None:
        raise ForbiddenHTTPException(f"Dataset {dataset_id} unlinked")
    if not user_info.is_user_owner_or_allowed(collection_version.owner):
        raise ForbiddenHTTPException("Unauthorized")
    # End of duplicate block

    # TODO: deduplicate from ApiCommon. We need to settle the class/module level debate before can do that
    if dataset_version.version_id not in [v.version_id for v in collection_version.datasets]:
        raise ForbiddenHTTPException(f"Dataset {dataset_id} does not belong to a collection")

    try:
        business_logic.remove_dataset_version(collection_version.version_id, dataset_version.version_id)
    except CollectionUpdateException:
        raise MethodNotAllowedException(detail="Cannot delete a public Dataset")
    return Response(status=202)
    # End of duplicate block


def put(collection_id: str, dataset_id: str, body: dict, token_info: dict):
    # TODO: deduplicate from ApiCommon. We need to settle the class/module level debate before can do that
    url = body.get("url", body.get("link"))
    business_logic = get_business_logic()
    dataset_version = business_logic.get_dataset_version(DatasetVersionId(dataset_id))
    if dataset_version is None:
        raise ForbiddenHTTPException(f"Dataset {dataset_id} does not exist")

    collection_version = business_logic.get_collection_version(CollectionVersionId(collection_id))
    if collection_version is None or not UserInfo(token_info).is_user_owner_or_allowed(collection_version.owner):
        raise ForbiddenHTTPException()

    try:
        business_logic.ingest_dataset(
            collection_version.version_id,
            url,
            None if dataset_id is None else DatasetVersionId(dataset_id),
        )
        return Response(status=202)
    except CollectionNotFoundException:
        raise ForbiddenHTTPException()
    except CollectionIsPublishedException:
        raise ForbiddenHTTPException()
    except DatasetNotFoundException:
        raise NotFoundHTTPException()
    except InvalidURIException:
        raise InvalidParametersHTTPException(detail="The dropbox shared link is invalid.")
    except MaxFileSizeExceededException:
        raise TooLargeHTTPException()
    except DatasetInWrongStatusException:
        raise MethodNotAllowedException(
            detail="Submission failed. A dataset cannot be updated while a previous update for the same dataset "
            "is in progress. Please cancel the current submission by deleting the dataset, or wait until "
            "the submission has finished processing."
        )
    # End of duplicate block
