import json
import logging
import os
from os import path
from typing import Dict, Optional

import dask
import numpy as np
import psutil
import tiledb
from cellxgene_schema.utils import read_h5ad

from backend.common.utils.corpora_constants import CorporaConstants
from backend.common.utils.cxg_constants import CxgConstants
from backend.common.utils.tiledb import consolidation_buffer_size
from backend.layers.processing.utils.color_conversion_utils import (
    ColorFormatException,
    convert_anndata_category_colors_to_cxg_category_colors,
)
from backend.layers.processing.utils.cxg_generation_utils import (
    convert_coverage_to_cxg_array,
    convert_dataframe_to_cxg_array,
    convert_dictionary_to_cxg_group,
    convert_matrices_to_cxg_arrays,
    convert_ndarray_to_cxg_dense_array,
    convert_uns_to_cxg_group,
)
from backend.layers.processing.utils.matrix_utils import is_matrix_sparse


class H5ADDataFile:
    """
    Class encapsulating required information about an H5AD datafile that ultimately will be transformed into
    another format (currently just CXG is supported).
    """

    tile_db_ctx_config = {
        "sm.consolidation.buffer_size": consolidation_buffer_size(0.1),
        "sm.consolidation.step_min_frags": 2,
        "sm.consolidation.step_max_frags": 20,  # see https://docs.tiledb.com/main/how-to/performance/performance-tips/tuning-consolidation
        "py.deduplicate": True,  # May reduce memory requirements at cost of performance
    }

    def __init__(
        self,
        input_filename,
        var_index_column_name=None,
    ):
        self.input_filename = input_filename
        # Set by self.extract_metadata_about_dataset
        self.dataset_title = None
        # These two are set by self.transform_dataframe_index_into_column
        self.var_index_column_name = var_index_column_name
        self.obs_index_column_name = None

        self.validate_input_file_type()

        self.extract_anndata_elements_from_file()
        self.extract_metadata_about_dataset()

        self.validate_anndata()

    def to_cxg(
        self,
        output_cxg_directory,
        sparse_threshold,
        dataset_version_id,
        fragment_artifact_id=None,
        convert_anndata_colors_to_cxg_colors=True,
    ):
        """
        Writes the following attributes of the anndata to CXG: 1) the metadata as metadata attached to an empty
        DenseArray, 2) the obs DataFrame as a DenseArray, 3) the var DataFrame as a DenseArray, 4) all valid
        embeddings stored in obsm, each one as a DenseArray, 5) the main X matrix of the anndata as either a
        SparseArray or DenseArray based on the `sparse_threshold`, and optionally 6) the column shift of the main X
        matrix that might turn an otherwise Dense matrix into a Sparse matrix.
        """

        logging.info("Beginning writing to CXG.")
        ctx = tiledb.Ctx(self.tile_db_ctx_config)

        tiledb.group_create(output_cxg_directory, ctx=ctx)
        logging.info(f"\t...group created, with name {output_cxg_directory}")

        convert_dictionary_to_cxg_group(
            output_cxg_directory, self.generate_cxg_metadata(convert_anndata_colors_to_cxg_colors), ctx=ctx
        )
        logging.info("\t...dataset metadata saved")

        convert_dataframe_to_cxg_array(output_cxg_directory, "obs", self.obs, self.obs_index_column_name, ctx)
        logging.info("\t...dataset obs dataframe saved")

        convert_dataframe_to_cxg_array(output_cxg_directory, "var", self.var, self.var_index_column_name, ctx)
        logging.info("\t...dataset var dataframe saved")

        convert_uns_to_cxg_group(output_cxg_directory, self.anndata.uns, dataset_version_id, "uns", ctx)
        logging.info("\t...dataset uns dataframe saved")

        if fragment_artifact_id is not None:
            deployment_stage = os.getenv("DEPLOYMENT_STAGE")
            if deployment_stage == "test":
                logging.info("Skipping ATAC processing in test environment")
            else:
                # Set the index of the obs dataframe to the obs_index_column_name
                # This is necessary for the coverage conversion to work correctly.
                # The obs_index_column_name is set in the transform_dataframe_index_into_column method.
                self.obs = self.obs.set_index(self.obs_index_column_name)

                convert_coverage_to_cxg_array(
                    output_cxg_directory, self.obs, fragment_artifact_id, "coverage", ctx, uns=self.anndata.uns
                )
                logging.info("\t...dataset coverage dataframe saved")

        self.write_anndata_embeddings_to_cxg(output_cxg_directory, ctx)
        logging.info("\t...dataset embeddings saved")

        self.write_anndata_x_matrices_to_cxg(output_cxg_directory, ctx, sparse_threshold)  # big memory usage
        logging.info("\t...dataset X matrix saved")

        logging.info("Completed writing to CXG.")

    def write_anndata_x_matrices_to_cxg(self, output_cxg_directory, ctx, sparse_threshold):
        matrix_container = f"{output_cxg_directory}/X"
        x_matrix_data = self.anndata.X
        with dask.config.set(
            {
                "num_workers": 1,  # a single worker with as many threads as vCPUs is more memory efficient
                "threads_per_worker": 2,  # match the number of vCPUs
                "distributed.worker.memory.limit": "14GB",
                "scheduler": "threads",
            }
        ):
            is_sparse = is_matrix_sparse(x_matrix_data, sparse_threshold)
            logging.info(f"is_sparse: {is_sparse}")

            convert_matrices_to_cxg_arrays(matrix_container, x_matrix_data, is_sparse, self.tile_db_ctx_config)

        logging.info("start consolidating")
        self._consolidate_tiledb_with_memory_optimization(matrix_container, ctx)

        if hasattr(tiledb, "vacuum"):
            tiledb.vacuum(matrix_container)

    def _consolidate_tiledb_with_memory_optimization(self, matrix_container: str, ctx: tiledb.Ctx):
        """Consolidate TileDB array with memory optimization."""
        available_memory_mb = psutil.virtual_memory().available / 1024 / 1024
        safe_buffer_mb = max(128, min(1024, int(available_memory_mb * 0.05)))
        safe_buffer_bytes = safe_buffer_mb * 1024 * 1024

        if available_memory_mb < 8192:
            logging.info(f"Low memory detected ({available_memory_mb:.0f}MB), using chunked consolidation approach")
            self._consolidate_tiledb_chunked(matrix_container, ctx, safe_buffer_bytes)
        else:
            self._consolidate_tiledb_standard(matrix_container, ctx, safe_buffer_bytes, safe_buffer_mb)

    def _consolidate_tiledb_standard(
        self, matrix_container: str, ctx: tiledb.Ctx, safe_buffer_bytes: int, safe_buffer_mb: int
    ):
        """Standard consolidation with memory optimization."""
        config_dict = {
            "sm.consolidation.buffer_size": safe_buffer_bytes,
            "sm.consolidation.steps_size_ratio": 0.1,
            "sm.memory_budget": safe_buffer_bytes * 2,
            "sm.consolidation.step_min_frags": 2,
            "sm.consolidation.step_max_frags": 10,
        }

        available_memory_mb = psutil.virtual_memory().available / 1024 / 1024
        logging.info(
            f"TileDB consolidation config: buffer={safe_buffer_mb}MB, "
            f"memory_budget={safe_buffer_mb * 2}MB, available_memory={available_memory_mb:.0f}MB"
        )

        try:
            config = tiledb.Config(config_dict)
            tiledb.consolidate(matrix_container, ctx=ctx, config=config)
            logging.info("TileDB consolidation completed successfully with memory optimization")
        except Exception as e:
            logging.warning(f"Memory-optimized consolidation failed: {e}")
            fallback_config_dict = {
                "sm.consolidation.buffer_size": safe_buffer_bytes // 2,
                "sm.memory_budget": safe_buffer_bytes,
                "sm.consolidation.steps_size_ratio": 0.05,
            }
            fallback_config = tiledb.Config(fallback_config_dict)
            logging.info(f"Retrying consolidation with fallback config: buffer={safe_buffer_mb//2}MB")
            tiledb.consolidate(matrix_container, ctx=ctx, config=fallback_config)

    def _consolidate_tiledb_chunked(self, matrix_container: str, ctx: tiledb.Ctx, safe_buffer_bytes: int):
        """Chunked consolidation approach that processes TileDB array subsets separately."""
        safe_buffer_mb = safe_buffer_bytes / (1024 * 1024)

        chunk_config_dict = {
            "sm.consolidation.buffer_size": safe_buffer_bytes // 4,
            "sm.memory_budget": safe_buffer_bytes // 2,
            "sm.consolidation.steps_size_ratio": 0.05,
            "sm.consolidation.step_min_frags": 2,
            "sm.consolidation.step_max_frags": 5,
        }

        logging.info(
            f"Chunked consolidation config: buffer={safe_buffer_mb/4:.0f}MB, " f"memory_budget={safe_buffer_mb/2:.0f}MB"
        )

        try:
            # Get array info to determine if chunking is beneficial
            array_info = tiledb.array_fragments(matrix_container, ctx=ctx)
            num_fragments = len(array_info)

            logging.info(f"Array has {num_fragments} fragments, consolidating in small batches")

            if num_fragments <= 5:
                # Few fragments, use standard approach with conservative settings
                chunk_config = tiledb.Config(chunk_config_dict)
                tiledb.consolidate(matrix_container, ctx=ctx, config=chunk_config)
                logging.info("Chunked consolidation completed (single pass)")
            else:
                # Many fragments, consolidate in multiple passes with fragment limits
                max_fragments_per_pass = 8
                passes_needed = max(1, (num_fragments + max_fragments_per_pass - 1) // max_fragments_per_pass)

                logging.info(
                    f"Performing {passes_needed} consolidation passes with max {max_fragments_per_pass} fragments per pass"
                )

                for pass_idx in range(passes_needed):
                    # Configure each pass to limit fragment count
                    pass_config_dict = chunk_config_dict.copy()
                    pass_config_dict["sm.consolidation.step_max_frags"] = max_fragments_per_pass
                    pass_config = tiledb.Config(pass_config_dict)

                    logging.info(f"Consolidation pass {pass_idx + 1}/{passes_needed}")
                    tiledb.consolidate(matrix_container, ctx=ctx, config=pass_config)

                    # Brief pause between passes to allow memory cleanup
                    import time

                    time.sleep(0.5)

                logging.info("Chunked consolidation completed (multi-pass)")

        except Exception as e:
            logging.warning(f"Chunked consolidation failed: {e}")
            # Last resort: minimal consolidation
            minimal_config_dict = {
                "sm.consolidation.buffer_size": safe_buffer_bytes // 8,
                "sm.memory_budget": safe_buffer_bytes // 4,
                "sm.consolidation.steps_size_ratio": 0.02,
                "sm.consolidation.step_max_frags": 3,
            }
            minimal_config = tiledb.Config(minimal_config_dict)
            logging.info("Attempting minimal consolidation as last resort")
            tiledb.consolidate(matrix_container, ctx=ctx, config=minimal_config)

    def write_anndata_embeddings_to_cxg(self, output_cxg_directory, ctx):
        def is_valid_embedding(adata, embedding_name, embedding_array):
            """
            Returns true if this layout data is a valid array for front-end presentation with the following criteria:
                * ndarray, with shape (n_obs, >= 2), dtype float/int/uint
                * follows ScanPy embedding naming conventions
                * with all values finite or NaN (no +Inf or -Inf)
            """

            is_valid = (
                isinstance(embedding_name, str)
                and (embedding_name.startswith("X_") or embedding_name == "spatial")
                and len(embedding_name) > 2
                and embedding_name != "X_spatial"
            )
            is_valid = is_valid and isinstance(embedding_array, np.ndarray) and embedding_array.dtype.kind in "fiu"
            is_valid = is_valid and embedding_array.shape[0] == adata.n_obs and embedding_array.shape[1] >= 2
            is_valid = is_valid and not np.any(np.isinf(embedding_array)) and not np.all(np.isnan(embedding_array))
            return is_valid

        embedding_container = f"{output_cxg_directory}/emb"
        tiledb.group_create(embedding_container, ctx=ctx)

        for embedding_name, embedding_values in self.anndata.obsm.items():
            if is_valid_embedding(self.anndata, embedding_name, embedding_values):
                if embedding_name == "spatial":
                    embedding_name = f"{embedding_container}/{embedding_name}"
                else:
                    embedding_name = f"{embedding_container}/{embedding_name[2:]}"
                convert_ndarray_to_cxg_dense_array(embedding_name, embedding_values, ctx)
                logging.info(f"\t\t...{embedding_name} embedding created")

    def generate_cxg_metadata(self, convert_anndata_colors_to_cxg_colors):
        """
        Return a dictionary containing metadata about CXG dataset. This include data about the version as well as
        Corpora schema properties if they exist, among other pieces of metadata.
        """

        cxg_group_metadata = {
            "cxg_version": CxgConstants.CXG_VERSION,
            "cxg_properties": json.dumps({"title": self.dataset_title, "about": None}),
        }
        if self.corpora_properties is not None:
            cxg_group_metadata["corpora"] = json.dumps(self.corpora_properties)

        if convert_anndata_colors_to_cxg_colors:
            try:
                cxg_group_metadata["cxg_category_colors"] = json.dumps(
                    convert_anndata_category_colors_to_cxg_category_colors(self.anndata)
                )
            except ColorFormatException:
                logging.warning(
                    "Failed to extract colors from H5AD file! Fix the H5AD file or rerun with "
                    "--disable-custom-colors. See help for more details."
                )

        return cxg_group_metadata

    def validate_input_file_type(self):
        """
        Validate that the input file is of a type that we can handle. Currently the only valid file type is `.h5ad`.
        """

        if not self.input_filename.endswith(".h5ad"):
            raise Exception(f"Cannot process input file {self.input_filename}. File must be an H5AD.")

    def validate_anndata(self):
        if not self.var.index.is_unique:
            raise ValueError("Variable index in AnnData object is not unique.")
        if not self.obs.index.is_unique:
            raise ValueError("Observation index in AnnData object is not unique.")

    def extract_anndata_elements_from_file(self):
        logging.info(f"Reading in AnnData dataset: {path.basename(self.input_filename)}")
        self.anndata = read_h5ad(self.input_filename, chunk_size=15000)
        logging.info("Completed reading in AnnData dataset!")

        self.obs = self.transform_dataframe_index_into_column(self.anndata.obs, "obs", self.obs_index_column_name)
        self.var = self.transform_dataframe_index_into_column(self.anndata.var, "var", self.var_index_column_name)

    def extract_metadata_about_dataset(self):
        """
        Extract metadata information about the dataset that upon conversion will be saved as group metadata with the
        CXG that is generated. This metadata information includes Corpora schema properties, the dataset title and
        a link that details more information about the dataset.
        """
        self.corpora_properties = self.get_corpora_properties()

        if self.corpora_properties is None:
            # If the return value is None, this means that we were not able to figure out what version of the Corpora
            # schema the object is using and therefore cannot extract any properties.
            raise ValueError("Unknown source file schema version is unsupported.")

        # The title and about properties of the dataset are set by the following order: if they are explicitly defined
        # then use the explicit value. Otherwise, use the input filename (only for title, about will be blank).

        filename = path.splitext(path.basename(self.input_filename))[0]
        self.dataset_title = self.dataset_title or filename

    def transform_dataframe_index_into_column(self, dataframe, dataframe_name, index_column_name):
        """
        Convert the dataframe's index into another column in the dataframe and reset the index to the default.
        If an index_column_name is specified, use that column as the index instead.
        """

        def convert_index_to_column():
            # Create a unique column name for the index.
            suffix = 0
            while f"name_{suffix}" in dataframe.columns:
                suffix += 1
            index_name = f"name_{suffix}"
            dataframe.rename_axis(index_name, inplace=True)
            dataframe.reset_index(inplace=True)
            return index_name

        if index_column_name is None:
            index_column_name = convert_index_to_column()

        elif index_column_name in dataframe.columns:
            # User has specified alternative column for unique names, and it exists
            if not dataframe[index_column_name].is_unique:
                raise KeyError(
                    f"Values in {dataframe_name}.{index_column_name} must be unique. Please prepare data to contain "
                    f"unique values."
                )
            # if not already in columns, save the current h5ad index as a column
            if dataframe.index.name not in dataframe.columns:
                convert_index_to_column()
        else:
            raise KeyError(f"Column {index_column_name} does not exist.")

        setattr(self, f"{dataframe_name}_index_column_name", index_column_name)
        return dataframe

    def get_corpora_properties(self) -> Optional[Dict]:
        """
        Extract out the Corpora dataset properties from the H5AD file.
        """
        corpora_props = {}
        for key in CorporaConstants.REQUIRED_SIMPLE_METADATA_FIELDS:
            if key not in self.anndata.uns:
                raise KeyError(f"missing Corpora schema field {key}")
            corpora_props[key] = self.anndata.uns[key]

        for key in CorporaConstants.OPTIONAL_LIST_METADATA_FIELDS:
            if key not in self.anndata.uns:
                continue
            try:
                corpora_props[key] = self.anndata.uns[key].tolist()
            except AttributeError as e:
                logging.error(f"Corpora schema field {key} expected to be list, got {type(self.anndata.uns[key])}")
                raise e

        for key in CorporaConstants.OPTIONAL_SIMPLE_METADATA_FIELDS:
            if key in self.anndata.uns:
                corpora_props[key] = self.anndata.uns[key]

        return corpora_props
