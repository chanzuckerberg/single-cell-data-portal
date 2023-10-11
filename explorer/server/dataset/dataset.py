from abc import ABCMeta, abstractmethod
from os.path import basename, splitext

import numpy as np
import pandas as pd
from server_timing import Timing as ServerTiming

from server.common.config.app_config import AppConfig
from server.common.constants import Axis, XApproximateDistribution
from server.common.errors import (
    DatasetAccessError,
    ExceedsLimitError,
    FilterError,
    JSONEncodingValueError,
    UnsupportedSummaryMethod,
)
from server.common.fbs.matrix import encode_matrix_fbs
from server.common.utils.utils import jsonify_numpy


class Dataset(metaclass=ABCMeta):
    """Base class for loading and accessing matrix data"""

    def __init__(self, data_locator, app_config):
        if not isinstance(app_config, AppConfig):
            raise TypeError("config expected to be of type AppConfig")

        # location to the dataset
        self.data_locator = data_locator

        # config is the application configuration
        self.app_config = app_config

        # parameters set by this data adaptor based on the data.
        self.parameters = {}

    @staticmethod
    @abstractmethod
    def pre_load_validation(data_locator):
        pass

    @staticmethod
    @abstractmethod
    def open(data_locator, app_config):
        pass

    @staticmethod
    @abstractmethod
    def file_size(data_locator):
        pass

    @abstractmethod
    def get_name(self):
        """return a string name for this data adaptor"""
        pass

    @abstractmethod
    def get_library_versions(self):
        """return a dictionary of library name to library versions"""
        pass

    @abstractmethod
    def get_embedding_names(self):
        """return a list of pre-computed embedding names"""
        pass

    @abstractmethod
    def get_embedding_array(self, ename, dims=2):
        """return an numpy array for the given pre-computed embedding name."""
        pass

    @abstractmethod
    def get_X_array(self, obs_mask=None, var_mask=None):
        """return the X array, possibly filtered by obs_mask or var_mask.
        the return type is ndarray."""
        pass

    @abstractmethod
    def get_X_approximate_distribution(self) -> XApproximateDistribution:
        """return the approximate distribution of the X matrix."""
        pass

    @abstractmethod
    def get_shape(self):
        pass

    @abstractmethod
    def query_var_array(self, term_var):
        pass

    @abstractmethod
    def query_obs_array(self, term_var):
        pass

    @abstractmethod
    def get_colors(self):
        pass

    @abstractmethod
    def get_obs_index(self):
        pass

    @abstractmethod
    def get_obs_columns(self):
        pass

    @abstractmethod
    def get_obs_keys(self):
        # return list of keys
        pass

    @abstractmethod
    def get_var_keys(self):
        # return list of keys
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cleanup()

    @abstractmethod
    def cleanup(self):
        pass

    def get_data_locator(self):
        return self.data_locator

    def get_location(self):
        return self.data_locator.uri_or_path

    def get_about(self):
        return None

    def get_title(self):
        # default to file name
        location = self.get_location()
        if location.endswith("/"):
            location = location[:-1]
        return splitext(basename(location))[0]

    def get_corpora_props(self):
        return None

    @abstractmethod
    def get_schema(self):
        """
        Return current schema
        """
        pass

    @abstractmethod
    def annotation_to_fbs_matrix(self, axis, field=None, uid=None, num_bins=None):
        """
        Gets annotation value for each observation
        :param axis: string obs or var
        :param fields: list of keys for annotation to return, returns all annotation values if not set.
        :param num_bins: number of bins for lossy integer compression. if None, no compression is performed.
        :return: flatbuffer: in fbs/matrix.fbs encoding
        """
        pass

    def update_parameters(self, parameters):
        parameters.update(self.parameters)

    def _index_filter_to_mask(self, filter, count):
        mask = np.zeros((count,), dtype=np.bool)
        for i in filter:
            if isinstance(i, list):
                mask[i[0] : i[1]] = True
            else:
                mask[i] = True
        return mask

    def _axis_filter_to_mask(self, axis, filter, count):
        mask = np.ones((count,), dtype=np.bool)
        if "index" in filter:
            mask = np.logical_and(mask, self._index_filter_to_mask(filter["index"], count))
        if "annotation_value" in filter:
            mask = np.logical_and(mask, self._annotation_filter_to_mask(axis, filter["annotation_value"], count))

        return mask

    def _annotation_filter_to_mask(self, axis, filter, count):
        mask = np.ones((count,), dtype=np.bool)
        for v in filter:
            name = v["name"]
            if axis == Axis.VAR:
                anno_data = self.query_var_array(name)
            elif axis == Axis.OBS:
                anno_data = self.query_obs_array(name)

            if anno_data.dtype.name in ["boolean", "category", "object"]:
                values = v.get("values", [])
                key_idx = np.in1d(anno_data, values)
                mask = np.logical_and(mask, key_idx)

            else:
                min_ = v.get("min", None)
                max_ = v.get("max", None)
                if min_ is not None:
                    key_idx = (anno_data >= min_).ravel()
                    mask = np.logical_and(mask, key_idx)
                if max_ is not None:
                    key_idx = (anno_data <= max_).ravel()
                    mask = np.logical_and(mask, key_idx)

        return mask

    def _filter_to_mask(self, filter):
        """
        Return the filter as a row and column selection list.
        No filter on a dimension means 'all'
        """
        shape = self.get_shape()
        var_selector = None
        obs_selector = None
        if filter is not None:
            if Axis.OBS in filter:
                obs_selector = self._axis_filter_to_mask(Axis.OBS, filter["obs"], shape[0])

            if Axis.VAR in filter:
                var_selector = self._axis_filter_to_mask(Axis.VAR, filter["var"], shape[1])

        return (obs_selector, var_selector)

    def data_frame_to_fbs_matrix(self, filter, axis, num_bins=None):
        """
        Retrieves data 'X' and returns in a flatbuffer Matrix.
        :param filter: filter: dictionary with filter params
        :param axis: string obs or var
        :param num_bins: number of bins for lossy integer compression. if None, no compression is performed.
        :return: flatbuffer Matrix

        Caveats:
        * currently only supports access on VAR axis
        * currently only supports filtering on VAR axis
        """
        with ServerTiming.time("where.query"):
            if axis != Axis.VAR:
                raise ValueError("Only VAR dimension access is supported")

            try:
                obs_selector, var_selector = self._filter_to_mask(filter)
            except (KeyError, IndexError, TypeError, AttributeError, DatasetAccessError):
                raise FilterError("Error parsing filter") from None

            if obs_selector is not None:
                raise FilterError("filtering on obs unsupported")

            num_columns = self.get_shape()[1] if var_selector is None else np.count_nonzero(var_selector)
            if self.app_config.exceeds_limit("column_request_max", num_columns):
                raise ExceedsLimitError("Requested dataframe columns exceed column request limit")

            X = self.get_X_array(obs_selector, var_selector)
        with ServerTiming.time("where.encode"):
            col_idx = np.nonzero([] if var_selector is None else var_selector)[0]
            fbs = encode_matrix_fbs(X, col_idx=col_idx, row_idx=None, num_bins=num_bins)

        return fbs

    def diffexp_topN(self, obsFilterA, obsFilterB, top_n=None):
        """
        Computes the top N differentially expressed variables between two observation sets. If mode
        is "TOP_N", then stats for the top N
        dataframes
        contain a subset of variables, then statistics for all variables will be returned, otherwise
        only the top N vars will be returned.
        :param obsFilterA: filter: dictionary with filter params for first set of observations
        :param obsFilterB: filter: dictionary with filter params for second set of observations
        :param top_n: Limit results to top N (Top var mode only)
        :return: top N genes and corresponding stats
        """
        if Axis.VAR in obsFilterA or Axis.VAR in obsFilterB:
            raise FilterError("Observation filters may not contain variable conditions")
        try:
            shape = self.get_shape()
            obs_mask_A = self._axis_filter_to_mask(Axis.OBS, obsFilterA["obs"], shape[0])
            obs_mask_B = self._axis_filter_to_mask(Axis.OBS, obsFilterB["obs"], shape[0])
        except (KeyError, IndexError):
            raise FilterError("Error parsing filter") from None
        if top_n is None:
            top_n = self.app_config.default_dataset__diffexp__top_n

        if self.app_config.exceeds_limit(
            "diffexp_cellcount_max", np.count_nonzero(obs_mask_A) + np.count_nonzero(obs_mask_B)
        ):
            raise ExceedsLimitError("Diffexp request exceeds max cell count limit")

        result = self.compute_diffexp_ttest(
            obs_mask_A,
            obs_mask_B,
            top_n=top_n,
            lfc_cutoff=self.app_config.default_dataset__diffexp__lfc_cutoff,
            selector_lists=False,
        )

        try:
            return jsonify_numpy(result)
        except ValueError:
            raise JSONEncodingValueError("Error encoding differential expression to JSON") from None

    def diffexp_topN_from_list(self, listA: np.ndarray, listB: np.ndarray, top_n: int = None):
        """
        Compute differential expression - same as diffexp_topN() - but specifying the
        two cell sets as lists of obs indices (postings lists).
        """
        if top_n is None:
            top_n = self.app_config.default_dataset__diffexp__top_n

        result = self.compute_diffexp_ttest(
            listA,
            listB,
            top_n=top_n,
            lfc_cutoff=self.app_config.default_dataset__diffexp__lfc_cutoff,
            selector_lists=True,
        )

        try:
            return jsonify_numpy(result)
        except ValueError:
            raise JSONEncodingValueError("Error encoding differential expression to JSON") from None

    @abstractmethod
    def compute_diffexp_ttest(self, maskA, maskB, top_n, lfc_cutoff, selector_lists=False):
        pass

    @staticmethod
    def normalize_embedding(embedding):
        """Normalize embedding layout to meet client assumptions.
        Embedding is an ndarray, shape (n_obs, n)., where n is normally 2
        """

        # scale isotropically
        try:
            min = np.nanmin(embedding, axis=0)
            max = np.nanmax(embedding, axis=0)
        except RuntimeError:
            # indicates entire array was NaN, which should propagate
            min = np.NaN
            max = np.NaN

        scale = np.amax(max - min)
        normalized_layout = (embedding - min) / scale

        # translate to center on both axis
        translate = 0.5 - ((max - min) / scale / 2)
        normalized_layout = normalized_layout + translate

        normalized_layout = normalized_layout.astype(dtype=np.float32)
        return normalized_layout

    def layout_to_fbs_matrix(self, fields, num_bins=None):
        """
        :param num_bins: number of bins for lossy integer compression. if None, no compression is performed.
        return specified embeddings as a flatbuffer, using the cellxgene matrix fbs encoding.

        * returns only first two dimensions, with name {ename}_0 and {ename}_1,
          where {ename} is the embedding name.
        * client assumes each will be individually centered & scaled (isotropically)
          to a [0, 1] range.
        * does not support filtering

        """
        embeddings = self.get_embedding_names() if fields is None or len(fields) == 0 else fields
        layout_data = []
        with ServerTiming.time("layout.query"):
            for ename in embeddings:
                embedding = self.get_embedding_array(ename, 2)
                normalized_layout = Dataset.normalize_embedding(embedding)
                layout_data.append(pd.DataFrame(normalized_layout, columns=[f"{ename}_0", f"{ename}_1"]))

        with ServerTiming.time("layout.encode"):
            df = pd.concat(layout_data, axis=1, copy=False) if layout_data else pd.DataFrame()
            fbs = encode_matrix_fbs(df, col_idx=df.columns, row_idx=None, num_bins=num_bins)

        return fbs

    def get_last_mod_time(self):
        try:
            lastmod = self.get_data_locator().lastmodtime()
        except RuntimeError:
            lastmod = None
        return lastmod

    def summarize_var(self, method, filter, query_hash, num_bins=None):
        with ServerTiming.time("summarize.query"):
            if method != "mean":
                raise UnsupportedSummaryMethod("Unknown gene set summary method.")

            obs_selector, var_selector = self._filter_to_mask(filter)
            if obs_selector is not None:
                raise FilterError("filtering on obs unsupported")

            # if no filter, just return zeros.  We don't have a use case
            # for summarizing the entire X without a filter, and it would
            # potentially be quite compute / memory intensive.
            if var_selector is None or np.count_nonzero(var_selector) == 0:
                mean = np.zeros((self.get_shape()[0], 1), dtype=np.float32)
            else:
                X = self.get_X_array(obs_selector, var_selector)
                mean = X.mean(axis=1, keepdims=True)
        with ServerTiming.time("summarize.encode"):
            col_idx = pd.Index([query_hash])
            fbs = encode_matrix_fbs(mean, col_idx=col_idx, row_idx=None, num_bins=num_bins)
        return fbs
