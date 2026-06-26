"""SPIKE: Parquet + DuckDB backend behind the CensusCubeQuery seam.

Throwaway proof-of-concept to measure parity + latency vs the TileDB cube. Mirrors the public
query methods of CensusCubeQuery (expression_summary, expression_summary_default, cell_counts,
marker_genes) but reads Parquet via DuckDB instead of TileDB arrays. Diffexp cubes are out of scope.

The TileDB-specific bits this replaces (query.py:_query):
  - indexed-dim slice  cube.df[(['ENSG..'],['UBERON..'],[])]   -> WHERE col IN (...)
  - predicate string   cube.query(cond="cell_type in [...]")    -> WHERE col IN (...)
  - attr/dim projection (reused verbatim via CensusCubeQueryParams + the schema shim below)
"""

from types import SimpleNamespace

import duckdb
import numpy as np
import pandas as pd

from backend.common.census_cube.data.query import depluralize


class _Col:
    def __init__(self, name, dtype):
        self.name = name
        self.dtype = np.dtype(dtype)


class _Schema:
    """Exposes just the surface CensusCubeQueryParams + _query introspect on a tiledb cube:
    iterate for attrs, .domain for indexed dims."""

    def __init__(self, domain, attrs):
        self.domain = domain
        self._attrs = attrs

    def __iter__(self):
        return iter(self._attrs)


class DuckDBCube:
    """A DuckDB table masquerading as a tiledb cube for the query layer."""

    def __init__(self, con: duckdb.DuckDBPyConnection, table: str, tiledb_schema):
        self.con = con
        self.table = table
        domain = [_Col(d.name, d.dtype) for d in tiledb_schema.domain]
        attrs = [_Col(a.name, a.dtype) for a in tiledb_schema]
        self.schema = _Schema(domain, attrs)
        self.columns = {c.name for c in domain} | {c.name for c in attrs}


def _sql_in(values):
    escaped = ", ".join("'" + str(v).replace("'", "''") + "'" for v in values)
    return f"({escaped})"


class DuckDBCensusCubeQuery:
    """Drop-in for CensusCubeQuery's read methods, backed by DuckDB/Parquet."""

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

    def _query(self, cube: DuckDBCube, criteria, compare_dimension=None) -> pd.DataFrame:
        # WHERE: every criteria field that maps to a column becomes an IN filter
        # (unifies tiledb's indexed-dim slice and non-indexed query-condition).
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

        select_cols = "*" if attrs is None and dims is None else ", ".join((dims or []) + (attrs or []))
        where = (" WHERE " + " AND ".join(clauses)) if clauses else ""
        return cube.con.execute(f"SELECT {select_cols} FROM {cube.table}{where}").df()


# --- store builder: fixture CSVs -> Parquet -> DuckDB ---------------------------------------


# (csv cube name, snapshot attr, tiledb schema) — schemas drive the dim/attr/dtype shim.
def build_duckdb_snapshot(snapshot_name: str, parquet_dir: str, fixtures_root: str) -> SimpleNamespace:
    from backend.common.census_cube.data.schemas.cube_schema import (
        cell_counts_schema,
        expression_summary_schema,
    )
    from backend.common.census_cube.data.schemas.cube_schema_default import (
        expression_summary_schema as expression_summary_default_schema,
    )
    from backend.common.census_cube.data.schemas.marker_gene_cube_schema import marker_genes_schema

    cubes = [
        ("expression_summary", "expression_summary_cube", expression_summary_schema),
        ("expression_summary_default", "expression_summary_default_cube", expression_summary_default_schema),
        ("cell_counts", "cell_counts_cube", cell_counts_schema),
        ("marker_genes", "marker_genes_cube", marker_genes_schema),
    ]

    con = duckdb.connect()
    snapshot = SimpleNamespace()
    for csv_name, attr, schema in cubes:
        df = pd.read_csv(f"{fixtures_root}/{snapshot_name}/{csv_name}.csv.gz", index_col=0)
        pq_path = f"{parquet_dir}/{csv_name}.parquet"
        df.to_parquet(pq_path, index=False)
        # Materialize once (in-memory table), the analog of TileDB's opened+cached array.
        # Swap to `CREATE VIEW ... read_parquet(...)` to measure cold scan-with-pushdown instead.
        con.execute(f"CREATE TABLE {csv_name} AS SELECT * FROM read_parquet('{pq_path}')")
        setattr(snapshot, attr, DuckDBCube(con, csv_name, schema))

    snapshot.cell_counts_df = pd.read_csv(f"{fixtures_root}/{snapshot_name}/cell_counts.csv.gz", index_col=0)
    snapshot._con = con
    return snapshot
