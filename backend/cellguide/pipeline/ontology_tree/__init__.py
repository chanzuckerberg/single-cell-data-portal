from backend.cellguide.pipeline.constants import (
    CELL_GUIDE_CELL_TYPE_MAPPING_FILENAME,
    ONTOLOGY_TREE_FILENAME,
    ONTOLOGY_TREE_STATE_PER_CELLTYPE_FOLDERNAME,
    ONTOLOGY_TREE_STATE_PER_TISSUE_FOLDERNAME,
    ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME,
)
from backend.cellguide.pipeline.ontology_tree.tree_builder import OntologyTreeBuilder
from backend.cellguide.pipeline.ontology_tree.types import OntologyTreeData
from backend.cellguide.pipeline.utils import output_json, output_json_per_key
from backend.wmg.api.wmg_api_config import WMG_API_SNAPSHOT_SCHEMA_VERSION
from backend.wmg.data.snapshot import WmgSnapshot, load_snapshot


def run(*, output_directory):
    snapshot = load_snapshot(snapshot_schema_version=WMG_API_SNAPSHOT_SCHEMA_VERSION)
    ontology_tree_data = get_ontology_tree_data(snapshot=snapshot)

    for organism in ontology_tree_data:
        celltype_to_tissue_mapping = get_celltype_to_tissue_mapping(ontology_tree_data[organism].all_states_per_tissue)

        organism_path_name = organism.replace(":", "_")
        output_json(
            ontology_tree_data[organism].ontology_graph,
            f"{output_directory}/{ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME}/{organism_path_name}/{ONTOLOGY_TREE_FILENAME}",
        )
        output_json_per_key(
            ontology_tree_data[organism].all_states_per_cell_type,
            f"{output_directory}/{organism_path_name}/{ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME}/{ONTOLOGY_TREE_STATE_PER_CELLTYPE_FOLDERNAME}",
        )
        output_json_per_key(
            ontology_tree_data[organism].all_states_per_tissue,
            f"{output_directory}/{organism_path_name}/{ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME}/{ONTOLOGY_TREE_STATE_PER_TISSUE_FOLDERNAME}",
        )
        output_json(
            celltype_to_tissue_mapping,
            f"{output_directory}/{organism_path_name}/{ONTOLOGY_TREE_TOPLEVEL_FOLDERNAME}/{CELL_GUIDE_CELL_TYPE_MAPPING_FILENAME}",
        )


def get_celltype_to_tissue_mapping(all_states_per_tissue):
    celltype_to_tissue_mapping = {}
    for tissue in all_states_per_tissue:
        tissueCounts = all_states_per_tissue[tissue].tissueCounts
        for celltype in tissueCounts:
            if tissueCounts[celltype]["n_cells_rollup"] > 0:
                L = celltype_to_tissue_mapping.get(celltype, set())
                L.add(tissue)
                celltype_to_tissue_mapping[celltype] = L
    celltype_to_tissue_mapping = {k: list(v) for k, v in celltype_to_tissue_mapping.items()}
    return celltype_to_tissue_mapping


def get_ontology_tree_builder(*, snapshot: WmgSnapshot) -> OntologyTreeBuilder:
    cell_counts_df = snapshot.cell_counts_cube.df[:]
    return OntologyTreeBuilder(cell_counts_df)


def get_ontology_tree_builder_for_organism(*, snapshot: WmgSnapshot, organism: str) -> OntologyTreeBuilder:
    cell_counts_df = snapshot.cell_counts_cube.df[:]
    cell_counts_df = cell_counts_df[cell_counts_df["organism_ontology_term_id"] == organism]
    return OntologyTreeBuilder(cell_counts_df)


def get_ontology_tree_data(*, snapshot: WmgSnapshot) -> tuple[OntologyTreeBuilder, OntologyTreeData]:
    organisms = snapshot.cell_counts_cube.df[:]["organism_ontology_term_id"].unique()
    ontology_tree_data = {}
    for organism in organisms:
        tree_builder = get_ontology_tree_builder_for_organism(snapshot=snapshot, organism=organism)
        ontology_graph = tree_builder.get_ontology_tree()
        all_states_per_cell_type = tree_builder.get_ontology_tree_state_per_celltype()
        all_states_per_tissue = tree_builder.get_ontology_tree_state_per_tissue()

        ontology_tree_data[organism] = OntologyTreeData(
            ontology_graph=ontology_graph,
            all_states_per_cell_type=all_states_per_cell_type,
            all_states_per_tissue=all_states_per_tissue,
        )

    return ontology_tree_data
