"""
Learning: TileDB push-down queries do not currently support conditions on variable length unicode strings.
It is a roadmap item.  Currently supported:
    * Unicode, fixed length
    * ASCII, variable length
Casting to ASCII for now as that covers 99.99% of our data (eg, ontology IDs).
"""

import pathlib
from collections import namedtuple
from typing import List, Union

import numpy as np
import pandas as pd
import tiledb

uint32_domain = (np.iinfo(np.uint32).min, np.iinfo(np.uint32).max - 1)

INTEGRATED_ARRAY_NAME = "integrated"
OBS_ARRAY_NAME = "obs"
VAR_ARRAY_NAME = "var"
DATASET_TO_GENE_IDS_NAME = "dataset_to_gene_ids"
FILTER_RELATIONSHIPS_NAME = "filter_relationships"


class LabelType(
    namedtuple(
        "Label",
        ["key", "dtype", "domain", "decode_from_index", "encode_as_dim", "var", "custom_decoder"],
        defaults=[None, False, False, None, None],
    )
):
    __slots__ = ()

    def _get_col(self, df: pd.DataFrame) -> Union[pd.Index, pd.Series, None]:
        if self.decode_from_index:
            return df.index
        elif self.key in df:
            return df[self.key]
        else:
            return None

    def decode(self, df: pd.DataFrame, *args, **kwargs) -> np.ndarray:
        if self.custom_decoder:
            return self.custom_decoder(self, df, *args, **kwargs)

        # else, simple dataframe column extraction
        col = self._get_col(df)
        if col is None:
            return np.full((df.shape[0],), b"", dtype="O")
        elif self.dtype is str:
            return col.to_numpy(dtype=str)
        elif self.dtype in ["ascii", np.bytes_]:
            return np.char.encode(col.to_numpy(dtype=str), encoding="ascii")
        else:
            return col.to_numpy(dtype=self.dtype)


def gen_idx(_lbl, df, start_coord=0):
    return np.arange(start_coord, start_coord + df.shape[0])


# dimensions - order matters for dimension
obs_labels = [
    # obs_idx is the join index with the X array, and IS ALSO the "iloc" from original dataset.
    LabelType("obs_idx", np.uint32, domain=uint32_domain, custom_decoder=gen_idx),
    LabelType("dataset_id", "ascii", encode_as_dim=True),
    LabelType("filter_cells", "bool", encode_as_dim=False),
    *[
        LabelType(key, "ascii", encode_as_dim=True)
        for key in [
            "cell_type_ontology_term_id",
            "tissue_ontology_term_id",
            "tissue_original_ontology_term_id",
        ]
    ],
    *[
        LabelType(key, "ascii", var=True)
        for key in [
            "cell_type",
            "assay",
            "assay_ontology_term_id",
            "development_stage",
            "development_stage_ontology_term_id",
            "disease_ontology_term_id",
            "tissue",
            "tissue_original",
            "self_reported_ethnicity",
            "self_reported_ethnicity_ontology_term_id",
            "sex",
            "sex_ontology_term_id",
            "organism",
            "organism_ontology_term_id",
        ]
    ],
    LabelType("dataset_local_cell_id", "ascii", var=True, decode_from_index=True),
]

# order matters for dimensions
var_labels = [
    # var_idx is the join index with the X array
    LabelType("var_idx", np.uint32, domain=uint32_domain, custom_decoder=gen_idx),
    LabelType("gene_ontology_term_id", "ascii", decode_from_index=True, encode_as_dim=True),
    LabelType("feature_reference", "ascii", var=True),
    LabelType("feature_name", "ascii", var=True),
]


def create_tdb_integrated_corpus(corpus_path: str):
    """
    Create the empty tiledb object for the integrated corpus
    """
    pathlib.Path(corpus_path).mkdir(parents=True, exist_ok=True)
    tiledb.group_create(corpus_path)

    filters = tiledb.FilterList([tiledb.ZstdFilter(level=-22)])
    create_integrated_expression_array(corpus_path, filters)
    create_obs_array(corpus_path, filters)
    create_var_array(corpus_path, filters)


def create_var_array(uri: str, filters: tiledb.filter.FilterList):
    """
    var/feature/gene axes labels.
    """
    tiledb.Array.create(
        f"{uri}/{VAR_ARRAY_NAME}",
        tiledb.ArraySchema(
            domain=tiledb.Domain(create_axes_label_dims(var_labels)),
            sparse=True,
            allows_duplicates=True,
            attrs=[
                tiledb.Attr(name=lbl.key, dtype=lbl.dtype, var=lbl.var, filters=filters)
                for lbl in var_labels
                if lbl.encode_as_dim is False
            ],
            cell_order="row-major",
            tile_order="row-major",
            capacity=10000,
        ),
    )


def create_obs_array(uri: str, filters: tiledb.filter.FilterList):
    """
    obs/cell axes labels
    """
    tiledb.Array.create(
        f"{uri}/{OBS_ARRAY_NAME}",
        tiledb.ArraySchema(
            domain=tiledb.Domain(create_axes_label_dims(obs_labels)),
            sparse=True,
            allows_duplicates=True,
            attrs=[
                tiledb.Attr(name=lbl.key, dtype=lbl.dtype, var=lbl.var, filters=filters)
                for lbl in obs_labels
                if lbl.encode_as_dim is False
            ],
            cell_order="row-major",
            tile_order="row-major",
            capacity=10000,
        ),
    )


def create_integrated_expression_array(uri: str, filters: tiledb.filter.FilterList):
    X_capacity = 128000
    X_extent = [512, 2048]  # guess - needs tuning
    tiledb.Array.create(
        f"{uri}/{INTEGRATED_ARRAY_NAME}",
        tiledb.ArraySchema(
            domain=tiledb.Domain(
                [
                    tiledb.Dim(
                        name="obs_idx",
                        domain=uint32_domain,
                        tile=X_extent[0],
                        dtype=np.uint32,
                        filters=filters,
                    ),
                    tiledb.Dim(
                        name="var_idx",
                        domain=uint32_domain,
                        tile=X_extent[1],
                        dtype=np.uint32,
                        filters=filters,
                    ),
                ]
            ),
            sparse=True,
            allows_duplicates=True,
            attrs=[
                tiledb.Attr(name="rankit", dtype=np.float32, filters=filters),
            ],
            cell_order="row-major",
            tile_order="col-major",
            capacity=X_capacity,
        ),
    )


def create_axes_label_dims(labels: List[LabelType]) -> List[tiledb.Dim]:
    dims = []
    extent = 1024  # guess - needs tuning
    filters = tiledb.FilterList([tiledb.ZstdFilter(level=22)])
    for lbl in labels:
        if not lbl.encode_as_dim:
            continue
        if lbl.encode_as_dim and lbl.dtype in [np.bytes_, "ascii"]:  # special case strings
            dim = tiledb.Dim(name=lbl.key, domain=None, tile=None, dtype=lbl.dtype, filters=filters, var=lbl.var)
        else:
            dim = tiledb.Dim(
                name=lbl.key, domain=lbl.domain, tile=extent, dtype=lbl.dtype, var=lbl.var, filters=filters
            )
        dims.append(dim)
    return dims
