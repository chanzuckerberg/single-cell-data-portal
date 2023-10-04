import logging
import math
import os
import pathlib
from typing import Optional, Tuple

import cellxgene_census
import numba as nb
import numpy as np
import pandas as pd
import tiledb
import tiledbsoma as soma
from cellxgene_census.experimental.util._csr_iter import X_sparse_iter

from backend.wmg.data.constants import (
    GENE_EXPRESSION_COUNT_MIN_THRESHOLD,
    RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD,
)
from backend.wmg.data.rankit import rankit
from backend.wmg.data.schemas.cube_schema import (
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
from backend.wmg.data.schemas.expression_summary_fmg_cube_schema import (
    expression_summary_fmg_indexed_dims,
    expression_summary_fmg_non_indexed_dims,
    expression_summary_fmg_schema,
)
from backend.wmg.data.snapshot import (
    EXPRESSION_SUMMARY_CUBE_NAME,
    EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME,
    EXPRESSION_SUMMARY_FMG_CUBE_NAME,
)
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.data.utils import (
    create_empty_cube,
    log_func_runtime,
    remove_accents,
    return_dataset_dict_w_publications,
)

logger = logging.getLogger(__name__)

DIMENSION_NAME_MAP_CENSUS_TO_WMG = {
    "tissue_ontology_term_id": "tissue_original_ontology_term_id",
    "tissue_general_ontology_term_id": "tissue_ontology_term_id",
}


@nb.njit(fastmath=True, error_model="numpy", parallel=False, nogil=True)
def gene_expression_sum_x_cube_dimension(
    rankit_values: np.ndarray,
    obs_idxs: np.ndarray,
    var_idx: np.ndarray,
    cube_indices: np.ndarray,
    sum_into: np.ndarray,
    sqsum_into: np.ndarray,
    nnz_into: np.ndarray,
):
    """
    Sum the rankit values for each gene (for each cube row/combo of cell attributes)
    Also track the number of cells that express that gene (nnz count)
    """
    for k in range(len(rankit_values)):
        val = rankit_values[k]
        if np.isfinite(val):
            cidx = var_idx[k]
            grp_idx = cube_indices[obs_idxs[k]]
            sum_into[grp_idx, cidx] += val
            sqsum_into[grp_idx, cidx] += val**2
            nnz_into[grp_idx, cidx] += 1


class SummaryCubesBuilder:
    def __init__(self, *, corpus_path: str, organismInfo: dict):
        """
        Initialize the SummaryCubesBuilder class.

        Args:
            corpus_path (str): The path to the corpus.
            organismInfo (dict): Information about the organism.
                ex:
                organismInfo = {
                    "id": "UBERON:0000955",
                    "label": "brain",
                }
        """
        organism = organismInfo["label"]
        organismId = organismInfo["id"]

        pathlib.Path(corpus_path).mkdir(exist_ok=True)
        with cellxgene_census.open_soma(census_version="latest") as census:
            obs_df = census["census_data"][organism]["obs"].read().concat().to_pandas().set_index("soma_joinid")
            obs_df["organism_ontology_term_id"] = organismId
            obs_df = obs_df[obs_df["is_primary_data"]]
            obs_df = obs_df.rename(columns=DIMENSION_NAME_MAP_CENSUS_TO_WMG)
            var_df = (
                census["census_data"][organism]["ms"]["RNA"]["var"].read().concat().to_pandas().set_index("soma_joinid")
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
        self.corpus_path = corpus_path
        self.organism = organism

        self.expression_summary_cube_created = False

    @log_func_runtime
    def create_expression_summary_cube(self, write_chunk_size: Optional[int] = 50_000_000):
        """
        Create the expression summary cube.

        Args:
            write_chunk_size (int, optional): The size of the write chunk. Defaults to 50_000_000.
        """
        uri = f"{self.corpus_path}/{EXPRESSION_SUMMARY_CUBE_NAME}"
        cube_dims = expression_summary_indexed_dims_no_gene_ontology + expression_summary_non_indexed_dims

        ctx = create_ctx()
        with tiledb.scope_ctx(ctx):
            if not os.path.exists(uri):
                create_empty_cube(uri, expression_summary_schema)

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

        self.expression_summary_cube_created = True

    @log_func_runtime
    def create_expression_summary_default_cube(self):
        """
        Create the default expression summary cube. The default expression summary cube is an aggregation across
        non-default dimensions in the expression summary cube.
        """
        logger.info("Creating the default expression summary cube.")
        if not self.expression_summary_cube_created:
            raise ValueError(
                "'expression_summary' array does not exist. Please run 'create_expression_summary_cube' first."
            )

        expression_summary_uri = f"{self.corpus_path}/{EXPRESSION_SUMMARY_CUBE_NAME}"
        expression_summary_default_uri = f"{self.corpus_path}/{EXPRESSION_SUMMARY_DEFAULT_CUBE_NAME}"

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

            if not os.path.exists(expression_summary_default_uri):
                create_empty_cube(expression_summary_default_uri, expression_summary_schema_default)
            logger.info(f"Writing cube to {expression_summary_default_uri}")
            tiledb.from_pandas(expression_summary_default_uri, expression_summary_df_default, mode="append")

    @log_func_runtime
    def create_expression_summary_fmg_cube(self):
        """
        Create the FMG (Find Marker Genes) expression summary cube. The FMG expression summary cube is an aggregation across
        dimensions in the expression summary cube that are irrelevant for marker gene calculation.
        """
        logger.info("Creating the FMG expression summary cube.")
        if not self.expression_summary_cube_created:
            raise ValueError(
                "'expression_summary' array does not exist. Please run 'create_expression_summary_cube' first."
            )

        expression_summary_uri = f"{self.corpus_path}/{EXPRESSION_SUMMARY_CUBE_NAME}"
        expression_summary_fmg_uri = f"{self.corpus_path}/{EXPRESSION_SUMMARY_FMG_CUBE_NAME}"

        ctx = create_ctx()
        with tiledb.scope_ctx(ctx):
            dfs = []
            with tiledb.open(expression_summary_uri, "r") as cube:
                for row in cube.query(return_incomplete=True).df[:]:
                    dfs.append(row)
            expression_summary_df = pd.concat(dfs, axis=0)

            expression_summary_df_fmg = (
                expression_summary_df.groupby(
                    expression_summary_fmg_indexed_dims + expression_summary_fmg_non_indexed_dims
                )
                .sum(numeric_only=True)
                .reset_index()
            )

            if not os.path.exists(expression_summary_fmg_uri):
                create_empty_cube(expression_summary_fmg_uri, expression_summary_fmg_schema)
            logger.info(f"Writing cube to {expression_summary_fmg_uri}")
            tiledb.from_pandas(expression_summary_fmg_uri, expression_summary_df_fmg, mode="append")

    @log_func_runtime
    def _summarize_gene_expressions(self, *, cube_dims: list, schema: tiledb.ArraySchema):
        """
        Summarize gene expressions for each row/combination of cell attributes.

        Args:
            cube_dims (list): The dimensions of the cube.
            schema (tiledb.ArraySchema): The schema of the cube.

        Returns:
            tuple: A tuple containing the dimensions and values of the cube.
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

            with organism.axis_query("RNA", obs_query=soma.AxisQuery(value_filter="is_primary_data == True")) as query:
                logger.info(f"Reducing X with {self.obs_df.shape[0]} total cells")

                iteration = 0
                for (obs_soma_joinids_chunk, _var_soma_joinids_chunk), raw_array in X_sparse_iter(
                    query, X_name="raw", stride=row_stride
                ):
                    logger.info(f"Reducer iteration {iteration} out of {math.ceil(self.obs_df.shape[0] / row_stride)}")
                    iteration += 1

                    obs_soma_joinids_chunk = query.indexer.by_obs(obs_soma_joinids_chunk)

                    gene_counts_per_cell = np.diff(raw_array.indptr)
                    keep = gene_counts_per_cell >= GENE_EXPRESSION_COUNT_MIN_THRESHOLD

                    # filter cells
                    raw_array = raw_array[keep]
                    obs_soma_joinids_chunk = obs_soma_joinids_chunk[keep]

                    data_filt = raw_array.data >= RANKIT_RAW_EXPR_COUNT_FILTERING_MIN_THRESHOLD

                    rankit(raw_array)

                    raw_array = raw_array.tocoo()

                    data = raw_array.data
                    row_indices, col = raw_array.row, raw_array.col
                    row = obs_soma_joinids_chunk[row_indices]

                    # filter data
                    keep_data = np.logical_and(np.isfinite(data), data_filt)
                    row, col, data = row[keep_data], col[keep_data], data[keep_data]

                    gene_expression_sum_x_cube_dimension(
                        rankit_values=data,
                        obs_idxs=row,
                        var_idx=col,
                        cube_indices=cube_indices,
                        sum_into=cube_sum,
                        nnz_into=cube_nnz,
                        sqsum_into=cube_sqsum,
                    )

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
