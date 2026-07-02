"""SPIKE: chDB (embedded ClickHouse MergeTree) backend behind the CensusCubeQuery seam.

Throwaway proof-of-concept to measure parity + latency vs the TileDB cube, mirroring query_lance.py.
The bet (vs Lance, which lost on a ~35ms row-id/bitmap floor): a MergeTree `ORDER BY (organism, gene,
tissue)` sparse primary index binary-searches granule marks for `WHERE gene IN (...)` and reads only
the matching 8192-row granules -> result-size access, structurally closer to TileDB's coordinate
tiling than Lance's scattered take. wmg_build_chdb.py builds the tables; the realcube benchmark
(scripts/wmg_chdb_spike_realcube.py) is what proves whether the granule pruning actually holds on the
351M-row cube -- same gate DuckDB and Lance failed.

Reuses the tiledb-schema shim (_Col/_Schema) and projection contract from query_duckdb so column
parity with the TileDB path holds by construction.
"""

import numpy as np
import pandas as pd

from backend.common.census_cube.data.query import depluralize
from backend.common.census_cube.data.query_duckdb import _Col, _Schema

DB = "wmg"


class ChdbCube:
    """A MergeTree table (in a chDB session) masquerading as a tiledb cube for the query layer."""

    def __init__(self, session, table: str, tiledb_schema):
        self.session = session
        self.table = table
        domain = [_Col(d.name, d.dtype) for d in tiledb_schema.domain]
        attrs = [_Col(a.name, a.dtype) for a in tiledb_schema]
        self.schema = _Schema(domain, attrs)
        self.columns = {c.name for c in domain} | {c.name for c in attrs}


def _sql_in(values):
    escaped = ", ".join("'" + str(v).replace("'", "''") + "'" for v in values)
    return f"({escaped})"


class ChdbCensusCubeQuery:
    """Drop-in for CensusCubeQuery's read methods, backed by chDB MergeTree tables."""

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

    def cell_counts_diffexp_df(self, criteria) -> pd.DataFrame:
        # Same pandas-filter as the TileDB path: cell_counts_diffexp is small and materialized in memory.
        df = self._snapshot.cell_counts_diffexp_df
        mask = np.ones(len(df), dtype=bool)
        for key, values in dict(criteria).items():
            values = values if isinstance(values, list) else [values]
            key = depluralize(key)
            if key in df.columns and values:
                mask &= df[key].isin(values)
        return df[mask].rename(columns={"n_cells": "n_total_cells"})

    def expression_summary_and_cell_counts_diffexp(self, criteria, use_simple: bool):
        # Mirror of CensusCubeQuery: get group_ids from the filtered cell_counts_diffexp, then slice the
        # diffexp cube by `group_id IN (...)`. MergeTree ORDER BY (group_id) binary-searches those groups.
        ccdf = self.cell_counts_diffexp_df(criteria)
        key = "group_id_simple" if use_simple else "group_id"
        cube = (
            self._snapshot.expression_summary_diffexp_simple_cube
            if use_simple
            else self._snapshot.expression_summary_diffexp_cube
        )
        group_ids = ccdf[key].unique().tolist()
        if not group_ids:
            return pd.DataFrame(columns=["group_id", "gene_ontology_term_id", "sum", "sqsum"]), ccdf
        in_list = "(" + ", ".join(str(int(g)) for g in group_ids) + ")"
        sql = (
            f"SELECT group_id, gene_ontology_term_id, sum, sqsum "
            f"FROM {DB}.{cube.table} WHERE group_id IN {in_list}"
        )
        return cube.session.query(sql, "DataFrame"), ccdf

    def _build_sql(self, cube: ChdbCube, criteria, compare_dimension=None):
        # WHERE: every criteria field that maps to a column becomes an IN filter. ClickHouse uses the
        # MergeTree primary key to binary-search granules for the leading dims (organism/gene) and
        # data-skipping indexes for the secondary categorical filters -> only matching granules read.
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

        columns = (dims or []) + (attrs or []) if (attrs is not None or dims is not None) else None
        select = ", ".join(columns) if columns else "*"
        where = " WHERE " + " AND ".join(clauses) if clauses else ""
        return f"SELECT {select} FROM {DB}.{cube.table}{where}"

    def _query(self, cube: ChdbCube, criteria, compare_dimension=None) -> pd.DataFrame:
        sql = self._build_sql(cube, criteria, compare_dimension)
        return cube.session.query(sql, "DataFrame")
