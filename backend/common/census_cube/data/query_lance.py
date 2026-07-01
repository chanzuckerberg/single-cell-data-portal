"""SPIKE: Lance backend behind the CensusCubeQuery seam.

Throwaway proof-of-concept to measure parity + latency vs the TileDB cube, mirroring query_duckdb.py.
The bet (vs DuckDB, which lost on full scans): Lance scalar indexes (BTREE on gene, BITMAP on the
low-cardinality dims) prune `col IN (...)` filters to result-size access instead of scanning — the
same layout property that makes TileDB's indexed-dim slicing fast. The realcube benchmark is what
proves whether that holds; see scripts/wmg_lance_spike_realcube.py.

Reuses the tiledb-schema shim (_Col/_Schema) and projection contract from query_duckdb so column
parity with the TileDB path holds by construction.
"""

import lance
import numpy as np
import pandas as pd

from backend.common.census_cube.data.query import depluralize
from backend.common.census_cube.data.query_duckdb import _Col, _Schema


class LanceCube:
    """A Lance dataset masquerading as a tiledb cube for the query layer."""

    def __init__(self, dataset: lance.LanceDataset, tiledb_schema):
        self.dataset = dataset
        domain = [_Col(d.name, d.dtype) for d in tiledb_schema.domain]
        attrs = [_Col(a.name, a.dtype) for a in tiledb_schema]
        self.schema = _Schema(domain, attrs)
        self.columns = {c.name for c in domain} | {c.name for c in attrs}


def _sql_in(values):
    escaped = ", ".join("'" + str(v).replace("'", "''") + "'" for v in values)
    return f"({escaped})"


class LanceCensusCubeQuery:
    """Drop-in for CensusCubeQuery's read methods, backed by Lance datasets."""

    def __init__(self, snapshot, cube_query_params=None):
        self._snapshot = snapshot
        self._cube_query_params = cube_query_params

    def expression_summary(self, criteria, compare_dimension=None):
        return self._query(self._snapshot.expression_summary_cube, criteria, compare_dimension)

    def expression_summary_default(self, criteria):
        return self._query(self._snapshot.expression_summary_default_cube, criteria)

    def marker_genes(self, criteria):
        return self._query(self._snapshot.marker_genes_cube, criteria)

    def cell_counts(self, criteria, compare_dimension=None):
        cc = self._query(
            self._snapshot.cell_counts_cube,
            criteria.copy(exclude={"gene_ontology_term_ids"}),
            compare_dimension,
        )
        return cc.rename(columns={"n_cells": "n_total_cells"})

    def _query(self, cube: LanceCube, criteria, compare_dimension=None) -> pd.DataFrame:
        # WHERE: every criteria field that maps to a column becomes an IN filter. Lance pushes these
        # into its scalar indexes (if built) -> only matching rows are materialized.
        clauses = []
        for field, vals in criteria.dict().items():
            col = depluralize(field)
            vals = vals if isinstance(vals, list) else [vals]
            vals = [v for v in vals if v != ""]
            if col in cube.columns and vals:
                clauses.append(f"{col} IN {_sql_in(vals)}")

        numeric_attrs = [a.name for a in cube.schema if np.issubdtype(a.dtype, np.number)]

        # Reuse the exact projection logic the tiledb path uses (column parity by construction).
        attrs = self._cube_query_params.get_attrs_for_cube_query(cube) if self._cube_query_params else None
        dims = self._cube_query_params.get_dims_for_cube_query(cube) if self._cube_query_params else None
        if attrs is not None:
            if compare_dimension is not None:
                attrs.append(compare_dimension)
            attrs += numeric_attrs

        columns = None if (attrs is None and dims is None) else (dims or []) + (attrs or [])
        where = " AND ".join(clauses) if clauses else None
        tbl = cube.dataset.to_table(columns=columns, filter=where)
        return tbl.to_pandas()
