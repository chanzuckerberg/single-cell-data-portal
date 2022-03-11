
"""
Learning: TileDB push-down queries do not currently support conditions on variable length unicode strings.
It is a roadmap item.  Currently supported:
    * Unicode, fixed length
    * ASCII, variable length
Casting to ASCII for now as that covers 99.99% of our data (eg, ontology IDs).
"""


# Hints on how to map between H5AD and TDB schemas.
from collections import namedtuple
from typing import Union

import numpy as np
import pandas as pd

uint32_domain = (np.iinfo(np.uint32).min, np.iinfo(np.uint32).max - 1)


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

    def decode(self, df, *args, **kwargs) -> np.ndarray:
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
    *[
        LabelType(key, "ascii", encode_as_dim=True)
        for key in [
            "cell_type_ontology_term_id",
            "tissue_ontology_term_id",
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
            "ethnicity",
            "ethnicity_ontology_term_id",
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


# order matters for dimensions
var_labels = [
    # var_idx is the join index with the X array
    LabelType("var_idx", np.uint32, domain=uint32_domain, custom_decoder=gen_idx),
    LabelType("gene_ontology_term_id", "ascii", decode_from_index=True, encode_as_dim=True),
    LabelType("feature_reference", "ascii", var=True),
    LabelType("feature_name", "ascii", var=True),
]