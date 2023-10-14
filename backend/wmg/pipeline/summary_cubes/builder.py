import json
import logging
import math
import os
import pathlib
from typing import Optional, Tuple

import cellxgene_census
import numpy as np
import pandas as pd
import tiledb
import tiledbsoma as soma
from cellxgene_census.experimental.util._csr_iter import X_sparse_iter

from backend.cellguide.pipeline.computational_marker_genes.computational_markers import MarkerGenesCalculator
from backend.cellguide.pipeline.ontology_tree import get_ontology_tree_builder
from backend.wmg.data.constants import (
    GENE_EXPRESSION_COUNT_MIN_THRESHOLD,
    RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD,
)
from backend.wmg.data.ontology_labels import gene_term_label, ontology_term_label
from backend.wmg.data.schemas.cube_schema import (
    cell_counts_logical_dims,
    cell_counts_schema,
    expression_summary_indexed_dims_no_gene_ontology,
    expression_summary_non_indexed_dims,
    expression_summary_schema,
)
from backend.wmg.data.schemas.cube_schema_default import (
    expression_summary_indexed_dims as expression_summary_indexed_dims_default,
)
from backend.wmg.data.schemas.cube_schema_default import (
    expression_summary_non_indexed_dims as expression_summary_non_indexed_dims_default,
)
from backend.wmg.data.schemas.cube_schema_default import expression_summary_schema as expression_summary_schema_default
from backend.wmg.data.schemas.marker_gene_cube_schema import marker_genes_schema
from backend.wmg.data.snapshot import (
    CELL_COUNTS_CUBE_NAME,
    DATASET_METADATA_FILENAME,
    EXPRESSION_SUMMARY_CUBE_NAME,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME,
    FILTER_RELATIONSHIPS_FILENAME,
    MARKER_GENES_CUBE_NAME,
    PRIMARY_FILTER_DIMENSIONS_FILENAME,
    WmgSnapshot,
)
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.utils import (
    create_empty_cube,
    get_datasets_from_curation_api,
    log_func_runtime,
)
from backend.wmg.pipeline.summary_cubes.builder_utils import (
    create_filter_relationships_graph,
    gene_expression_sum_x_cube_dimension,
    list_grouped_primary_filter_dimensions_term_ids,
    order_tissues,
    remove_accents,
    return_dataset_dict_w_publications,
)

logger = logging.getLogger(__name__)

DIMENSION_NAME_MAP_CENSUS_TO_WMG = {
    "tissue_ontology_term_id": "tissue_original_ontology_term_id",
    "tissue_general_ontology_term_id": "tissue_ontology_term_id",
}

EXPRESSION_SUMMARY_CUBE_CREATED_FLAG = "expression_summary_cube_created"
CELL_COUNTS_CUBE_CREATED_FLAG = "cell_counts_cube_created"
EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG = "expression_summary_default_cube_created"
MARKER_GENES_CUBE_CREATED_FLAG = "marker_genes_cube_created"
FILTER_RELATIONSHIPS_CREATED_FLAG = "filter_relationships_created"
PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG = "primary_filter_dimensions_created"
DATASET_METADATA_CREATED_FLAG = "dataset_metadata_created"


class SummaryCubesBuilder:
    def __init__(self, *, corpus_path: str, organismInfo: dict):
        """
        Initialize the SummaryCubesBuilder class.

        Args:
            corpus_path (str): The path to the corpus.
            organismInfo (dict): Information about the organism.
                ex:
                organismInfo = {
                    "id": "NCBITaxon:9606",
                    "label": "homo_sapiens",
                }
                Note that "label" should be the organism label used by census.
        """
        self.snapshot_id = os.path.basename(os.path.normpath(corpus_path))

        pathlib.Path(corpus_path).mkdir(parents=True, exist_ok=True)
        state_file = os.path.join(corpus_path, "pipeline_state.json")
        if os.path.exists(state_file):
            with open(state_file, "r") as f:
                self.state = json.load(f)
        else:
            self.state = {
                EXPRESSION_SUMMARY_CUBE_CREATED_FLAG: False,
                CELL_COUNTS_CUBE_CREATED_FLAG: False,
                EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG: False,
                MARKER_GENES_CUBE_CREATED_FLAG: False,
                FILTER_RELATIONSHIPS_CREATED_FLAG: False,
                PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG: False,
                DATASET_METADATA_CREATED_FLAG: False,
            }
        self.state_file = state_file

        self.organism = organismInfo["label"]
        self.organismId = organismInfo["id"]
        self.corpus_path = os.path.join(corpus_path, organismInfo["id"].replace(":", "_"))

        logger.info(f"Creating directory {self.corpus_path}")
        pathlib.Path(self.corpus_path).mkdir(parents=True, exist_ok=True)

    @log_func_runtime
    def _load_obs_and_var_dfs_if_necessary(self):
        """
        This method loads the obs and var dataframes from the census data.

        Raises:
            ValueError: If the organism is not a valid Census organism.
        """
        if hasattr(self, "obs_df") and hasattr(self, "var_df") and hasattr(self, "census_version"):
            return

        logger.info("Loading obs and var dataframes from the census data...")
        with cellxgene_census.open_soma(census_version="latest") as census:
            valid_organisms = list(census["census_data"].keys())
            if self.organism not in valid_organisms:
                raise ValueError(
                    f'"{self.organism}" is not a valid Census organism. Please select one of {valid_organisms}'
                )

            obs_df = census["census_data"][self.organism]["obs"].read().concat().to_pandas().set_index("soma_joinid")
            obs_df = obs_df[obs_df["is_primary_data"] & (obs_df["nnz"] >= GENE_EXPRESSION_COUNT_MIN_THRESHOLD)]
            obs_df = obs_df.rename(columns=DIMENSION_NAME_MAP_CENSUS_TO_WMG)
            var_df = (
                census["census_data"][self.organism]["ms"]["RNA"]["var"]
                .read()
                .concat()
                .to_pandas()
                .set_index("soma_joinid")
            )

            self.census_version = (
                census["census_info"]["summary"]
                .read()
                .concat()
                .to_pandas()
                .set_index("label")["value"]["census_build_date"]
            )
        self.obs_df = obs_df
        self.var_df = var_df

    @log_func_runtime
    def run_pipeline(self):
        """
        This method runs the entire pipeline for creating summary cubes. It sequentially calls the methods for creating
        the expression summary cube, the default expression summary cube, the cell counts cube and filter relationships,
        the primary filter dimensions, and the marker genes cube.
        """
        self.create_expression_summary_cube()
        self.create_expression_summary_default_cube()
        self.create_cell_counts_cube_and_filter_relationships()
        self.create_primary_filter_dimensions()
        self.create_marker_genes_cube()
        self.create_dataset_metadata()

    @log_func_runtime
    def create_expression_summary_default_cube(self):
        """
        Create the default expression summary cube. The default expression summary cube is an aggregation across
        non-default dimensions in the expression summary cube.
        """
        if not self._is_step_completed(EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG):
            raise ValueError(
                "'expression_summary' array does not exist. Please run 'create_expression_summary_cube' first."
            )

        logger.info("Creating the default expression summary cube.")
        expression_summary_uri = os.path.join(self.corpus_path, EXPRESSION_SUMMARY_CUBE_NAME)
        expression_summary_default_uri = os.path.join(self.corpus_path, EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME)

        ctx = create_ctx()
        with tiledb.scope_ctx(ctx):
            dfs = []
            with tiledb.open(expression_summary_uri, "r") as cube:
                for row in cube.query(return_incomplete=True).df[:]:
                    dfs.append(row)
            expression_summary_df = pd.concat(dfs, axis=0)

            expression_summary_df_default = (
                expression_summary_df.groupby(
                    expression_summary_indexed_dims_default + expression_summary_non_indexed_dims_default
                )
                .sum(numeric_only=True)
                .reset_index()
            )

            self._create_empty_cube(expression_summary_default_uri, expression_summary_schema_default)
            logger.info(f"Writing cube to {expression_summary_default_uri}")
            tiledb.from_pandas(expression_summary_default_uri, expression_summary_df_default, mode="append")

            self._mark_step_completed(EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG)

    @log_func_runtime
    def create_cell_counts_cube_and_filter_relationships(self):
        """
        Create cell count cube and write to disk
        """
        self._load_obs_and_var_dfs_if_necessary()

        logger.info("Creating the cell counts cube and filter relationships graph.")

        df = (
            self.obs_df.groupby(
                by=[dim for dim in cell_counts_logical_dims if dim != "publication_citation"],
                as_index=False,
            ).size()
        ).rename(columns={"size": "n_cells"})

        dataset_dict = return_dataset_dict_w_publications()
        df["publication_citation"] = [
            remove_accents(dataset_dict.get(dataset_id, "No Publication")) for dataset_id in df["dataset_id"]
        ]

        n_cells = df["n_cells"].to_numpy()
        df["n_cells"] = n_cells

        logger.info("Creating and writing filter relationships graph.")
        filter_relationships_linked_list = create_filter_relationships_graph(df)
        with open(f"{self.corpus_path}/{FILTER_RELATIONSHIPS_FILENAME}", "w") as f:
            json.dump(filter_relationships_linked_list, f)

        self._mark_step_completed(FILTER_RELATIONSHIPS_CREATED_FLAG)

        uri = os.path.join(self.corpus_path, CELL_COUNTS_CUBE_NAME)
        self._create_empty_cube(uri, cell_counts_schema)
        logger.info("Writing cell counts cube.")
        tiledb.from_pandas(uri, df, mode="append")

        self._mark_step_completed(CELL_COUNTS_CUBE_CREATED_FLAG)

    @log_func_runtime
    def create_marker_genes_cube(self):
        expression_summary_default_cube_uri = os.path.join(self.corpus_path, EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME)
        cell_counts_cube_uri = os.path.join(self.corpus_path, CELL_COUNTS_CUBE_NAME)
        if not self._is_step_completed(EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG):
            raise ValueError(
                "'expression_summary_default' cube does not exist. Please run 'create_expression_summary_default_cube' first."
            )

        if not self._is_step_completed(CELL_COUNTS_CUBE_CREATED_FLAG):
            raise ValueError(
                "'cell_counts' cube does not exist. Please run 'create_cell_counts_cube_and_filter_relationships' first."
            )

        if not self._is_step_completed(PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG):
            raise ValueError(
                "'primary_filter_dimensions' file does not exist. Please run 'create_primary_filter_dimensions' first."
            )

        logger.info("Calculating marker genes.")
        with (
            open(os.path.join(self.corpus_path, PRIMARY_FILTER_DIMENSIONS_FILENAME), "r") as f,
            tiledb.open(cell_counts_cube_uri, "r") as cell_counts_cube,
            tiledb.open(expression_summary_default_cube_uri, "r") as expression_summary_default_cube,
        ):
            primary_filter_dimensions = json.load(f)
            snapshot = WmgSnapshot(
                primary_filter_dimensions=primary_filter_dimensions,
                cell_counts_cube=cell_counts_cube,
                expression_summary_default_cube=expression_summary_default_cube,
            )
            ontology_tree = get_ontology_tree_builder(snapshot=snapshot)
            calculator = MarkerGenesCalculator(
                snapshot=snapshot,
                all_cell_type_ids_in_corpus=ontology_tree.all_cell_type_ids_in_corpus,
                groupby_terms=["tissue_ontology_term_id"],
            )
            marker_genes = calculator.get_computational_marker_genes()
        marker_gene_records = []
        for key in marker_genes:
            for marker_gene in marker_genes[key]:
                marker_gene_records.append(
                    {
                        "tissue_ontology_term_id": marker_gene.groupby_dims["tissue_ontology_term_id"],
                        "cell_type_ontology_term_id": key,
                        "gene_ontology_term_id": marker_gene.gene_ontology_term_id,
                        "marker_score": marker_gene.marker_score,
                        "specificity": marker_gene.specificity,
                    }
                )
        marker_genes_df = pd.DataFrame(marker_gene_records)
        uri = os.path.join(self.corpus_path, MARKER_GENES_CUBE_NAME)
        self._create_empty_cube(uri, marker_genes_schema)
        logger.info("Writing marker genes cube.")
        tiledb.from_pandas(uri, marker_genes_df, mode="append")

        self._mark_step_completed(MARKER_GENES_CUBE_CREATED_FLAG)

    @log_func_runtime
    def create_expression_summary_cube(self, write_chunk_size: Optional[int] = 50_000_000):
        """
        Create the expression summary cube.

        Args:
            write_chunk_size (int, optional): The size of the write chunk. Defaults to 50_000_000.
        """
        self._load_obs_and_var_dfs_if_necessary()

        uri = os.path.join(self.corpus_path, EXPRESSION_SUMMARY_CUBE_NAME)
        cube_dims = expression_summary_indexed_dims_no_gene_ontology + expression_summary_non_indexed_dims

        ctx = create_ctx()
        with tiledb.scope_ctx(ctx):
            self._create_empty_cube(uri, expression_summary_schema)
            dims, vals = self._summarize_gene_expressions(cube_dims=cube_dims, schema=expression_summary_schema)

            num_chunks = math.ceil(dims[0].size / write_chunk_size)

            with tiledb.open(uri, "w") as cube:
                for i in range(num_chunks):
                    logger.info(f"Writing chunk {i} to {uri}")
                    dims_chunk = tuple(dim[i * write_chunk_size : (i + 1) * write_chunk_size] for dim in dims)
                    vals_chunk = {k: vals[k][i * write_chunk_size : (i + 1) * write_chunk_size] for k in vals}
                    cube[dims_chunk] = vals_chunk

            logger.info("Consolidating and vacuuming")
            tiledb.consolidate(uri)
            tiledb.vacuum(uri)

            self._mark_step_completed(EXPRESSION_SUMMARY_CUBE_CREATED_FLAG)

    @log_func_runtime
    def create_primary_filter_dimensions(self):
        """
        This method creates the primary filter dimensions for the WMG snapshot.
        """
        if not self._is_step_completed(EXPRESSION_SUMMARY_DEFAULT_CUBE_CREATED_FLAG):
            raise ValueError(
                "'expression_summary_default' array does not exist. Please run 'create_expression_summary_default_cube' first."
            )
        if not self._is_step_completed(CELL_COUNTS_CUBE_CREATED_FLAG):
            raise ValueError(
                "'cell_counts' array does not exist. Please run 'create_cell_counts_cube_and_filter_relationships' first."
            )

        logger.info("Creating primary filter dimensions")
        with (
            tiledb.open(f"{self.corpus_path}/{EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME}", "r") as expression_summary_cube,
            tiledb.open(f"{self.corpus_path}/{CELL_COUNTS_CUBE_NAME}", "r") as cell_counts_cube,
        ):
            # get dataframes
            cell_counts_df = cell_counts_cube.df[:]
            expr_df = expression_summary_cube.df[:]

            # populate organism ID into the dfs for the primary filter dimensions
            cell_counts_df["organism_ontology_term_id"] = self.organismId
            expr_df["organism_ontology_term_id"] = self.organismId

            # genes
            organism_gene_ids = list_grouped_primary_filter_dimensions_term_ids(
                expr_df, "gene_ontology_term_id", "organism_ontology_term_id"
            )
            organism_gene_terms = {
                organism_term_id: [{g: gene_term_label(g)} for g in gene_term_ids]
                for organism_term_id, gene_term_ids in organism_gene_ids.items()
            }

            # tissues
            organism_tissue_ids = list_grouped_primary_filter_dimensions_term_ids(
                cell_counts_df, "tissue_ontology_term_id", group_by_dim="organism_ontology_term_id"
            )

            # organisms
            organism_tissue_terms = {
                organism_term_id: [{t: ontology_term_label(t)} for t in order_tissues(tissue_term_ids)]
                for organism_term_id, tissue_term_ids in organism_tissue_ids.items()
            }

            # collate
            result = dict(
                snapshot_id=str(self.snapshot_id),
                organism_terms=[
                    {o: ontology_term_label(o)} for o in sorted(cell_counts_df["organism_ontology_term_id"].unique())
                ],
                tissue_terms=organism_tissue_terms,
                gene_terms=organism_gene_terms,
            )
            logger.info("Writing primary filter dimensions")
            with open(f"{self.corpus_path}/{PRIMARY_FILTER_DIMENSIONS_FILENAME}", "w") as f:
                json.dump(result, f)

            self._mark_step_completed(PRIMARY_FILTER_DIMENSIONS_CREATED_FLAG)

    @log_func_runtime
    def create_dataset_metadata(self, snapshot_path: str) -> None:
        """
        This function generates a dictionary containing metadata for each dataset.
        The metadata includes the dataset id, label, collection id, and collection label.
        The function fetches the datasets from the curation API and iterates over them to create the metadata dictionary.
        """
        logger.info("Generating dataset metadata file")
        datasets = get_datasets_from_curation_api()
        dataset_dict = {}
        for dataset in datasets:
            dataset_id = dataset["dataset_id"]
            dataset_dict[dataset_id] = dict(
                id=dataset_id,
                label=dataset["title"],
                collection_id=dataset["collection_id"],
                collection_label=dataset["collection_name"],
            )
        logger.info("Writing dataset metadata file")
        with open(f"{snapshot_path}/{DATASET_METADATA_FILENAME}", "w") as f:
            json.dump(dataset_dict, f)
        self._mark_step_completed(DATASET_METADATA_CREATED_FLAG)

    @log_func_runtime
    def _summarize_gene_expressions(self, *, cube_dims: list, schema: tiledb.ArraySchema):
        """
        Summarize gene expressions for each row/combination of cell attributes.

        Args:
            cube_dims (list): The dimensions of the cube.
            schema (tiledb.ArraySchema): The schema of the cube.

        Returns:
            tuple: A tuple containing the dimensions (list) and values (keys) of the cube.
        """
        cube_index, cell_labels = self._make_cube_index(
            cube_dims=[dim for dim in cube_dims if dim != "publication_citation"]
        )

        n_groups = len(cube_index)
        n_genes = len(self.var_df)

        logger.info(f"Summarizing gene expressions across {n_groups} groups and {n_genes} genes")

        cube_sum = np.zeros((n_groups, n_genes), dtype=np.float32)
        cube_nnz = np.zeros((n_groups, n_genes), dtype=np.uint64)
        cube_sqsum = np.zeros((n_groups, n_genes), dtype=np.float32)

        self._reduce_X(
            cube_indices=cell_labels.cube_idx.values,
            cube_sum=cube_sum,
            cube_nnz=cube_nnz,
            cube_sqsum=cube_sqsum,
        )

        dim_names = [dim.name for dim in schema.domain]
        dims, vals = self._build_in_mem_cube(
            schema=schema,
            cube_index=cube_index,
            other_cube_attrs=[i for i in cube_dims if i not in dim_names],
            cube_sum=cube_sum,
            cube_nnz=cube_nnz,
            cube_sqsum=cube_sqsum,
        )
        return dims, vals

    @log_func_runtime
    def _make_cube_index(self, *, cube_dims: list) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Create the cube index. The cube index is a dataframe containing all possible groups of metadata with
        their corresponding indices.

        Args:
            cube_dims (list): The dimensions of the cube.

        Returns:
            Tuple[pd.DataFrame, pd.DataFrame]: A tuple containing the cube index and cell labels.
        """
        logger.info(f"Making the cube index across {cube_dims}")
        cell_labels = pd.DataFrame(
            data={k: self.obs_df[k].astype("category") for k in cube_dims},
            index=self.obs_df.index,
        )

        # number of cells with specific tuple of dims
        cube_index = pd.DataFrame(cell_labels.value_counts(), columns=["n"])

        # add cube_idx column
        cube_index["cube_idx"] = range(len(cube_index))
        cube_index["cube_idx"] = cube_index["cube_idx"].astype("int")

        # join cube_idx to cell_labels
        cell_labels = cell_labels.join(cube_index.cube_idx, on=cube_dims)
        cell_labels["cube_idx"] = cell_labels["cube_idx"].astype("int")

        return cube_index, cell_labels

    @log_func_runtime
    def _reduce_X(
        self,
        *,
        cube_indices: np.ndarray,
        cube_sum: np.ndarray,
        cube_nnz: np.ndarray,
        cube_sqsum: np.ndarray,
    ) -> None:
        """
        This method reduces the X matrix and stores the results in cube_sum, cube_nnz, and cube_sqsum.
        """
        with cellxgene_census.open_soma(
            census_version=self.census_version,
        ) as census:
            organism = census["census_data"][self.organism]
            row_stride = 50_000

            with organism.axis_query(
                "RNA",
                obs_query=soma.AxisQuery(
                    value_filter=f"is_primary_data == True and nnz >= {GENE_EXPRESSION_COUNT_MIN_THRESHOLD}"
                ),
            ) as query:
                logger.info(f"Reducing X with {self.obs_df.shape[0]} total cells")

                iteration = 0
                for (obs_soma_joinids_chunk, _), raw_array in X_sparse_iter(
                    query,
                    X_name="normalized",
                    stride=row_stride,
                    fmt="csr",
                ):
                    logger.info(f"Reducer iteration {iteration} out of {math.ceil(self.obs_df.shape[0] / row_stride)}")
                    iteration += 1

                    # convert soma_joinids to row coords of filtered obs dataframe
                    obs_soma_joinids_chunk = query.indexer.by_obs(obs_soma_joinids_chunk)

                    # convert to coo TODO: remove once X_sparse_iter can return coo
                    raw_array = raw_array.tocoo()

                    # get data and coordinates
                    data = raw_array.data
                    row_indices, col_var = raw_array.row, raw_array.col

                    # convert row coordinates to obs coordinates
                    row_obs = obs_soma_joinids_chunk[row_indices]

                    # convert data to raw counts and filter by min threshold
                    data_filt = (
                        self.obs_df["raw_sum"].values[row_obs] * raw_array.data
                        >= RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD
                    )

                    # convert data to log2(CPM+1)
                    raw_array.data[:] = np.log2(raw_array.data * 1e6 + 1)

                    # filter data
                    keep_data = np.logical_and(np.isfinite(data), data_filt)
                    row_obs, col_var, data = row_obs[keep_data], col_var[keep_data], data[keep_data]

                    gene_expression_sum_x_cube_dimension(
                        rankit_values=data,
                        obs_idxs=row_obs,
                        var_idx=col_var,
                        cube_indices=cube_indices,
                        sum_into=cube_sum,
                        nnz_into=cube_nnz,
                        sqsum_into=cube_sqsum,
                    )

    @log_func_runtime
    def _build_in_mem_cube(
        self,
        *,
        schema: tiledb.ArraySchema,
        cube_index: pd.DataFrame,
        other_cube_attrs: list,
        cube_sum: np.ndarray,
        cube_nnz: np.ndarray,
        cube_sqsum: np.ndarray,
    ):
        """
        Build the cube in memory, calculating the gene expression value for each combination of attributes
        """
        logger.info("Building cube in memory")
        # Count total values so we can allocate buffers once
        total_vals = 0
        for cube_idx in cube_index.cube_idx.values:
            mask = cube_nnz[cube_idx] != 0
            total_vals += np.count_nonzero(mask)

        # allocate buffers
        dims = [np.empty((total_vals,), dtype=object) for i in range(len(schema.domain))]

        vals = {
            "sum": np.empty((total_vals,)),
            "sqsum": np.empty((total_vals,)),
            "nnz": np.empty((total_vals,), dtype=np.uint64),
            **{k: np.empty((total_vals,), dtype=object) for k in other_cube_attrs},
        }

        # populate buffers
        idx = 0

        if "publication_citation" in other_cube_attrs:
            dataset_dict = return_dataset_dict_w_publications()

        for grp in cube_index.to_records():
            (
                *dim_or_attr_values,
                _,
                cube_idx,
            ) = grp.tolist()

            dim_or_attr_values = dict(zip(cube_index.index.names, dim_or_attr_values))

            mask = cube_nnz[cube_idx] != 0
            n_vals = np.count_nonzero(mask)
            if n_vals == 0:  # Used to maintain sparsity
                continue

            for i, dim in enumerate(schema.domain):
                if dim.name == "gene_ontology_term_id":
                    dims[i][idx : idx + n_vals] = self.var_df.feature_id.values[mask]
                else:
                    dims[i][idx : idx + n_vals] = dim_or_attr_values[dim.name]

            vals["sum"][idx : idx + n_vals] = cube_sum[cube_idx, mask]
            vals["nnz"][idx : idx + n_vals] = cube_nnz[cube_idx, mask]
            vals["sqsum"][idx : idx + n_vals] = cube_sqsum[cube_idx, mask]

            for _i, k in enumerate(other_cube_attrs):
                if k != "publication_citation":
                    vals[k][idx : idx + n_vals] = dim_or_attr_values[k]

            if "publication_citation" in other_cube_attrs:
                vals["publication_citation"][idx : idx + n_vals] = remove_accents(
                    dataset_dict.get(dim_or_attr_values["dataset_id"], "No Publication")
                )

            idx += n_vals

        return dims, vals

    def _create_empty_cube(self, uri: str, schema: tiledb.ArraySchema):
        logger.info(f"Creating empty cube at {uri}")
        create_empty_cube(uri, schema)

    def _mark_step_completed(self, step_key: str):
        """
        Mark the pipeline step complete for a given key and save the state to disk.

        Args:
            step_key (str): The key in the state dictionary to update.
        """
        self.state[step_key] = True
        with open(self.state_file, "w") as f:
            json.dump(self.state, f)

    def _is_step_completed(self, step_key: str) -> bool:
        """
        Check if a pipeline step is completed.

        Args:
            step_key (str): The key in the state dictionary to check.

        Returns:
            bool: True if the step is completed, False otherwise.
        """
        return self.state.get(step_key, False)
