from backend.cellguide.pipeline.canonical_marker_genes.utils import format_citation_dp
from backend.cellguide.pipeline.source_collections.types import SourceCollectionsData
from backend.common.utils.rollup import descendants
from backend.wmg.data.utils import get_collections_from_curation_api, get_datasets_from_curation_api


def generate_source_collections_data(all_cell_type_ids_in_corpus: list[str]) -> dict[str, list[SourceCollectionsData]]:
    """
    For each cell type id in the corpus, we want to generate a SourceCollectionsData object, which contains
    metadata about the source data for each cell type
    """
    all_datasets = get_datasets_from_curation_api()
    all_collections = get_collections_from_curation_api()

    collections_dict = {collection["collection_id"]: collection for collection in all_collections}
    datasets_dict = {dataset["dataset_id"]: dataset for dataset in all_datasets}

    source_collections_data: dict[str, list[SourceCollectionsData]] = {}

    for cell_id in all_cell_type_ids_in_corpus:
        lineage = descendants(cell_id)
        assert cell_id in lineage, f"{cell_id} not found in lineage"

        # Generate a list of unique dataset ids that contain cell type ids in the lineage
        dataset_ids: list[str] = []
        for dataset_id in datasets_dict:
            for cell_type in datasets_dict[dataset_id]["cell_type"]:
                if cell_type["ontology_term_id"] in lineage:
                    dataset_ids.append(dataset_id)
                    break
        assert len(set(dataset_ids)) == len(dataset_ids)

        # Generate a mapping from collection id to SourceCollectionsData
        collections_to_source_data = {}
        for dataset_id in dataset_ids:
            dataset = datasets_dict[dataset_id]
            collection_id = dataset["collection_id"]

            # If we don't have any source data on this collection yet, create a new SourceCollectionsData object
            if collection_id not in collections_to_source_data:
                collection = collections_dict.get(collection_id)
                source_data = SourceCollectionsData(
                    collection_name=collection["name"],
                    collection_url=collection["collection_url"],
                    publication_url=collection["doi"],
                    publication_title=format_citation_dp(collection["publisher_metadata"])
                    if collection["publisher_metadata"]
                    else "No publication",
                    tissue=[],
                    disease=[],
                    organism=[],
                )
                collections_to_source_data[collection_id] = source_data

            # Add the tissue, disease, and organism metadata from the dataset to the SourceCollectionsData object. If we
            # previously found a SourceCollectionsData object, we'll want to add the tissue/disease/organism from this
            # dataset to the existing object. If not, we'll want to create a new list for each of these fields.
            source_data = collections_to_source_data[collection_id]

            for tissue in dataset["tissue"]:
                if tissue["ontology_term_id"] not in [t["ontology_term_id"] for t in source_data.tissue]:
                    source_data.tissue.append(tissue)

            for disease in dataset["disease"]:
                if disease["ontology_term_id"] not in [d["ontology_term_id"] for d in source_data.disease]:
                    source_data.disease.append(disease)

            for organism in dataset["organism"]:
                if organism["ontology_term_id"] not in [o["ontology_term_id"] for o in source_data.organism]:
                    source_data.organism.append(organism)

            # Add the SourceCollectionsData to the mapping
            collections_to_source_data[collection_id] = source_data

        source_collections_data[cell_id] = list(collections_to_source_data.values())

    return source_collections_data
