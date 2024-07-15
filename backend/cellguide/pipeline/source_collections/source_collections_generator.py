import cellxgene_census
from pandas import DataFrame

from backend.cellguide.pipeline.source_collections.types import SourceCollectionsData
from backend.common.census_cube.utils import descendants


def generate_source_collections_data(
    all_cell_type_ids_in_corpus: list[str], cell_counts_df: DataFrame
) -> dict[str, list[SourceCollectionsData]]:
    """
    For each cell type id in the corpus, we want to generate a SourceCollectionsData object, which contains
    metadata about the source data for each cell type
    """
    dataset_id_to_cell_type_ids_map = {}
    dataset_id_to_tissue_ids_map = {}
    dataset_id_to_disease_ids_map = {}
    dataset_id_to_organism_ids_map = {}

    for map_dict, column_name in zip(
        [
            dataset_id_to_cell_type_ids_map,
            dataset_id_to_tissue_ids_map,
            dataset_id_to_disease_ids_map,
            dataset_id_to_organism_ids_map,
        ],
        [
            "cell_type_ontology_term_id",
            "tissue_ontology_term_id",
            "disease_ontology_term_id",
            "organism_ontology_term_id",
        ],
        strict=False,
    ):
        df_agg = cell_counts_df.groupby("dataset_id").agg({column_name: lambda x: ",".join(set(x.values))})
        df_dict = {df_agg.index[i]: df_agg.values[i][0].split(",") for i in range(len(df_agg))}
        map_dict.update(df_dict)

    with cellxgene_census.open_soma(census_version="latest") as census:
        datasets_metadata_df = census["census_info"]["datasets"].read().concat().to_pandas()

    source_collections_data: dict[str, list[SourceCollectionsData]] = {}
    for cell_id in all_cell_type_ids_in_corpus:
        lineage = descendants(cell_id)
        assert cell_id in lineage, f"{cell_id} not found in lineage"

        # Generate a list of unique dataset ids that contain cell type ids in the lineage
        dataset_ids: list[str] = []
        for dataset_id in dataset_id_to_cell_type_ids_map:
            for cell_type in dataset_id_to_cell_type_ids_map[dataset_id]:
                if cell_type in lineage:
                    dataset_ids.append(dataset_id)
                    break
        assert len(set(dataset_ids)) == len(dataset_ids)

        # Generate a mapping from collection id to SourceCollectionsData
        collections_to_source_data = {}

        for i in range(len(datasets_metadata_df)):
            dataset_id = datasets_metadata_df.iloc[i]["dataset_id"]
            if dataset_id not in dataset_ids:
                continue

            assert dataset_id in cell_counts_df["dataset_id"].unique(), f"{dataset_id} not in cell_counts_df"

            collection_id = datasets_metadata_df.iloc[i]["collection_id"]
            if collection_id not in collections_to_source_data:
                source_data = SourceCollectionsData(
                    collection_name=datasets_metadata_df.iloc[i]["collection_name"],
                    collection_url=f"https://cellxgene.cziscience.com/collections/{collection_id}",
                    publication_url=datasets_metadata_df.iloc[i]["collection_doi"],
                    publication_title=datasets_metadata_df.iloc[i]["collection_doi_label"],
                    tissue=dataset_id_to_tissue_ids_map[dataset_id],
                    disease=dataset_id_to_disease_ids_map[dataset_id],
                    organism=dataset_id_to_organism_ids_map[dataset_id],
                )
                collections_to_source_data[collection_id] = source_data

        source_collections_data[cell_id] = list(collections_to_source_data.values())
    return source_collections_data
