import logging
import math
import os

import numpy as np
import pandas as pd
import tiledb
from numba import njit
from scipy import sparse
from tiledbsoma import ExperimentAxisQuery

from backend.wmg.data.schemas.cube_schema import (
    expression_summary_indexed_dims_no_gene_ontology,
    expression_summary_non_indexed_dims,
    expression_summary_schema,
)
from backend.wmg.data.snapshot import EXPRESSION_SUMMARY_CUBE_NAME
from backend.wmg.data.tiledb import create_ctx
from backend.wmg.pipeline.constants import (
    ASSAYS_FOR_GENE_LENGTH_NORMALIZATION,
    DIMENSION_NAME_MAP_CENSUS_TO_WMG,
    NORM_EXPR_COUNT_FILTERING_MIN_THRESHOLD,
    TARGET_LIBRARY_SIZE,
)
from backend.wmg.pipeline.utils import (
    create_empty_cube_if_needed,
    load_pipeline_state,
    log_func_runtime,
    remove_accents,
    return_dataset_dict_w_publications,
)

logger = logging.getLogger(__name__)


WRITE_CHUNK_SIZE = 50_000_000


class ExpressionSummaryCubeBuilder:
    def __init__(self, *, query: ExperimentAxisQuery, corpus_path: str, organismId: str):
        self.obs_df = query.obs().concat().to_pandas()
        self.obs_df = self.obs_df.rename(columns=DIMENSION_NAME_MAP_CENSUS_TO_WMG)
        self.obs_df["organism_ontology_term_id"] = organismId

        self.var_df = query.var().concat().to_pandas()
        self.query = query
        self.corpus_path = corpus_path

        self.pipeline_state = load_pipeline_state(corpus_path)

    @log_func_runtime
    def create_expression_summary_cube(self):
        uri = os.path.join(self.corpus_path, EXPRESSION_SUMMARY_CUBE_NAME)
        cube_dims = expression_summary_indexed_dims_no_gene_ontology + expression_summary_non_indexed_dims

        ctx = create_ctx()
        with tiledb.scope_ctx(ctx):
            create_empty_cube_if_needed(uri, expression_summary_schema)
            dims, vals = self._summarize_gene_expressions(cube_dims=cube_dims, schema=expression_summary_schema)

            num_chunks = math.ceil(dims[0].size / WRITE_CHUNK_SIZE)

            with tiledb.open(uri, "w") as cube:
                for i in range(num_chunks):
                    logger.info(f"Writing chunk {i} to {uri}")
                    dims_chunk = tuple(dim[i * WRITE_CHUNK_SIZE : (i + 1) * WRITE_CHUNK_SIZE] for dim in dims)
                    vals_chunk = {k: vals[k][i * WRITE_CHUNK_SIZE : (i + 1) * WRITE_CHUNK_SIZE] for k in vals}
                    cube[dims_chunk] = vals_chunk

            logger.info("Consolidating and vacuuming")
            tiledb.consolidate(uri)
            tiledb.vacuum(uri)

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
    def _make_cube_index(self, *, cube_dims: list) -> tuple[pd.DataFrame, pd.DataFrame]:
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
        row_stride = 50_000

        logger.info(f"Reducing X with {self.obs_df.shape[0]} total cells")

        iteration = 0
        for raw_array, (obs_soma_joinids_chunk, _) in self.query.X("raw").blockwise(axis=0, size=row_stride).scipy():
            assert isinstance(raw_array, sparse.csr_matrix)
            logger.info(f"Reducer iteration {iteration} out of {math.ceil(self.obs_df.shape[0] / row_stride)}")
            iteration += 1

            # convert soma_joinids to row coords of filtered obs dataframe
            obs_soma_joinids_chunk = self.query.indexer.by_obs(obs_soma_joinids_chunk)

            # get assay ids for gene length normalization
            assay_ids = self.obs_df["assay_ontology_term_id"].values[obs_soma_joinids_chunk]
            # get boolean array of which cells to normalize
            cells_for_gene_length_normalization = np.in1d(assay_ids, ASSAYS_FOR_GENE_LENGTH_NORMALIZATION)

            num_cells_for_norm = cells_for_gene_length_normalization.sum()

            # note that `.multiply` converts the array to COO format
            if num_cells_for_norm == len(obs_soma_joinids_chunk):
                # normalize by gene length
                raw_array = raw_array.multiply(1 / self.var_df["feature_length"].values[None, :])
                # normalize by library size
                raw_array = raw_array.multiply(TARGET_LIBRARY_SIZE / raw_array.sum(1).A)

            elif num_cells_for_norm == 0:
                # normalize by library size
                raw_array = raw_array.multiply(TARGET_LIBRARY_SIZE / raw_array.sum(1).A)
            else:
                # get indices of cells to not normalize
                no_norm_idx = np.where(~cells_for_gene_length_normalization)[0]
                # normalized subset of raw_array by library size
                raw_arr_no_norm = raw_array[no_norm_idx]
                raw_arr_no_norm = raw_arr_no_norm.multiply(TARGET_LIBRARY_SIZE / raw_arr_no_norm.sum(1).A)

                # get indices of cells to normalize
                with_norm_idx = np.where(cells_for_gene_length_normalization)[0]
                # normalized subset of raw_array by gene length and library size
                raw_arr_with_norm = raw_array[with_norm_idx]
                raw_arr_with_norm = raw_arr_with_norm.multiply(1 / self.var_df["feature_length"].values[None, :])
                raw_arr_with_norm = raw_arr_with_norm.multiply(TARGET_LIBRARY_SIZE / raw_arr_with_norm.sum(1).A)

                # get indexers to remap row indices to their corresponding row indices in raw_array
                no_norm_indexer = pd.Series(index=np.arange(raw_arr_no_norm.shape[0]), data=no_norm_idx)
                with_norm_indexer = pd.Series(index=np.arange(raw_arr_with_norm.shape[0]), data=with_norm_idx)

                # get the row and column indices of the normalized and non-normalized arrays
                no_norm_row_idx, no_norm_col_idx = raw_arr_no_norm.row, raw_arr_no_norm.col
                with_norm_row_idx, with_norm_col_idx = raw_arr_with_norm.row, raw_arr_with_norm.col

                # remap row indices to their corresponding row indices in raw_array
                no_norm_row_idx = no_norm_indexer[no_norm_row_idx].values
                with_norm_row_idx = with_norm_indexer[with_norm_row_idx].values

                # combine row and column indices and data from the two arrays
                combined_row_idx = np.append(no_norm_row_idx, with_norm_row_idx)
                combined_col_idx = np.append(no_norm_col_idx, with_norm_col_idx)
                combined_data = np.append(raw_arr_no_norm.data, raw_arr_with_norm.data)

                # recreate raw_array
                raw_array = sparse.coo_matrix(
                    (combined_data, (combined_row_idx, combined_col_idx)), shape=raw_array.shape
                )

            # get data and coordinates
            data = raw_array.data
            row_indices, col_var = raw_array.row, raw_array.col

            # convert row coordinates to obs coordinates
            row_obs = obs_soma_joinids_chunk[row_indices]

            # convert data to raw counts and filter by min threshold
            data_filt = data >= NORM_EXPR_COUNT_FILTERING_MIN_THRESHOLD

            # convert data to log(X+1)
            data = np.log(data + 1)

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

            for _, k in enumerate(other_cube_attrs):
                if k != "publication_citation":
                    vals[k][idx : idx + n_vals] = dim_or_attr_values[k]

            if "publication_citation" in other_cube_attrs:
                vals["publication_citation"][idx : idx + n_vals] = remove_accents(
                    dataset_dict.get(dim_or_attr_values["dataset_id"], "No Publication")
                )

            idx += n_vals

        return dims, vals


@njit(fastmath=True, error_model="numpy", parallel=False, nogil=True)
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
