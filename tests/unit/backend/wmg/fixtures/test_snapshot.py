import contextlib
import json
import os
import sys
import tempfile
import uuid
from collections import namedtuple
from itertools import cycle, filterfalse, islice
from typing import Callable, Dict, List, NamedTuple, Tuple

import numpy as np
import pandas as pd
import tiledb
from numpy.random import randint, random
from pandas import DataFrame

from backend.wmg.data.schemas.cube_schema import (
    cell_counts_indexed_dims,
    cell_counts_logical_attrs,
    cell_counts_logical_dims,
    cell_counts_schema,
    expression_summary_indexed_dims,
    expression_summary_logical_attrs,
    expression_summary_logical_dims,
    expression_summary_schema,
)
from backend.wmg.data.snapshot import CELL_TYPE_ORDERINGS_FILENAME, WmgSnapshot
from backend.wmg.data.tiledb import create_ctx
from tests.unit.backend.wmg.fixtures import FIXTURES_ROOT
from tests.unit.backend.wmg.fixtures.test_primary_filters import build_precomputed_primary_filters


def simple_ontology_terms_generator(dimension_name: str, n_terms: int) -> List[str]:
    return [f"{dimension_name}_{i}" for i in range(n_terms)]


def semi_real_dimension_values_generator(dimension_name: str, dim_size: int) -> List[str]:
    """
    Returns a set of ontology term ids, sampled from real ontologies. While these ontology terms are
    from the real ontologies, they are not necessarily ones that would be admissible w.r.t. to the cellxgene dataset
    schema. This is still useful to ensure that id-to-label mapping lookups will return a real label. Note that this
    implementation is wildly inefficient, but it is good enough for test code.
    """
    # must import lazily
    import backend.wmg.data.ontology_labels as ontology_labels

    if ontology_labels.ontology_term_id_labels is None:
        ontology_labels.__load_ontologies()
    if ontology_labels.gene_term_id_labels is None:
        ontology_labels.__load_genes()

    deterministic_term_ids = sorted(ontology_labels.ontology_term_id_labels.keys())

    if dimension_name == "gene_ontology_term_id":
        return sorted(ontology_labels.gene_term_id_labels.keys())[:dim_size]
    if dimension_name == "tissue_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("UBERON")][:dim_size]
    if dimension_name == "tissue_original_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("UBERON")][:dim_size]
    if dimension_name == "organism_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("NCBITaxon")][:dim_size]
    if dimension_name == "cell_type_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("CL")][:dim_size]
    if dimension_name == "dataset_id":
        return [str(uuid.UUID()) for i in range(dim_size)]
    if dimension_name == "assay_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("EFO")][:dim_size]
    if dimension_name == "development_stage_ontology_term_id":
        return [
            term_id for term_id in deterministic_term_ids if term_id.startswith("Hsap") or term_id.startswith("MmusDev")
        ][:dim_size]
    if dimension_name == "disease_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("MONDO")][:dim_size]
    if dimension_name == "self_reported_ethnicity_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("HANCESTRO")][:dim_size]
    if dimension_name == "sex_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("PATO")][:dim_size]
    raise AssertionError(f"unknown dimension name {dimension_name}")


def random_expression_summary_values(coords):
    return {
        "nnz": randint(size=len(coords), low=0, high=100),
        "sum": random(size=len(coords)) * 10,
    }


def all_ones_expression_summary_values(coords):
    return {
        "nnz": np.ones(len(coords)),
        "sum": np.ones(len(coords)),
    }


def all_tens_cell_counts_values(coords) -> List[int]:
    return list(np.full(shape=len(coords), fill_value=10.0))


def all_X_cell_counts_values(coords, X) -> List[int]:
    return list(np.full(shape=len(coords), fill_value=X))


def random_cell_counts_values(coords) -> List[int]:
    return list(randint(size=len(coords), low=1, high=1000))


def exclude_random_coords_75pct(_) -> bool:
    return random() > 0.75


def exclude_dev_stage_and_ethnicity_for_secondary_filter_test(coord) -> bool:
    dev_stages_to_exclude = ("development_stage_ontology_term_id_1", "development_stage_ontology_term_id_2")
    self_reported_ethnicity_terms_to_exclude = (
        "self_reported_ethnicity_ontology_term_id_1",
        "self_reported_ethnicity_ontology_term_id_2",
    )
    if (
        coord.development_stage_ontology_term_id in dev_stages_to_exclude
        and coord.self_reported_ethnicity_ontology_term_id in self_reported_ethnicity_terms_to_exclude
    ):
        return True
    return False


# use this to create a disjoint set of genes across organisms, to mimic real data (each organism has its own set of
# genes); without this filtering function, the cube would have the cross-product of organisms * genes
# noinspection PyUnresolvedReferences
def exclude_all_but_one_gene_per_organism(logical_coord: NamedTuple) -> bool:
    # HACK: method called during building of both "expr summary" and "cell count" cubes, but the latter does not
    # include gene_ontology_term_id
    if "gene_ontology_term_id" not in logical_coord._fields:
        return False
    return logical_coord.gene_ontology_term_id != logical_coord.organism_ontology_term_id.replace("organism", "gene")


def forward_cell_type_ordering(cell_type_ontology_ids: List[str]) -> List[int]:
    return list(range(len(cell_type_ontology_ids)))


def reverse_cell_type_ordering(cell_type_ontology_ids: List[str]) -> List[int]:
    return list(range(len(cell_type_ontology_ids), 0, -1))


@contextlib.contextmanager
def load_realistic_test_snapshot(snapshot_name: str) -> WmgSnapshot:
    with tiledb.open(
        f"{FIXTURES_ROOT}/{snapshot_name}/expression_summary", ctx=create_ctx()
    ) as expression_summary_cube, tiledb.open(
        f"{FIXTURES_ROOT}/{snapshot_name}/expression_summary_default", ctx=create_ctx()
    ) as expression_summary_default_cube, tiledb.open(
        f"{FIXTURES_ROOT}/{snapshot_name}/expression_summary_fmg", ctx=create_ctx()
    ) as expression_summary_fmg_cube, tiledb.open(
        f"{FIXTURES_ROOT}/{snapshot_name}/cell_counts", ctx=create_ctx()
    ) as cell_counts_cube, tiledb.open(
        f"{FIXTURES_ROOT}/{snapshot_name}/marker_genes", ctx=create_ctx()
    ) as marker_genes_cube, open(
        f"{FIXTURES_ROOT}/{snapshot_name}/dataset_to_gene_ids.json", "r"
    ) as f:
        dataset_to_gene_ids = json.load(f)
        yield WmgSnapshot(
            snapshot_identifier=snapshot_name,
            expression_summary_cube=expression_summary_cube,
            expression_summary_fmg_cube=expression_summary_fmg_cube,
            expression_summary_default_cube=expression_summary_default_cube,
            marker_genes_cube=marker_genes_cube,
            cell_counts_cube=cell_counts_cube,
            cell_type_orderings=None,
            primary_filter_dimensions=None,
            dataset_to_gene_ids=dataset_to_gene_ids,
        )


@contextlib.contextmanager
def create_temp_wmg_snapshot(
    dim_size=3,
    snapshot_name="dummy-snapshot",
    expression_summary_vals_fn: Callable[[List[Tuple]], Dict[str, List]] = random_expression_summary_values,
    exclude_logical_coord_fn: Callable[[NamedTuple], bool] = None,
    cell_counts_generator_fn: Callable[[List[Tuple]], List] = random_cell_counts_values,
    cell_ordering_generator_fn: Callable[[List[str]], List[int]] = forward_cell_type_ordering,
) -> WmgSnapshot:
    with tempfile.TemporaryDirectory() as cube_dir:
        expression_summary_cube_dir, cell_counts_cube_dir = create_cubes(
            cube_dir,
            dim_size,
            exclude_logical_coord_fn=exclude_logical_coord_fn,
            expression_summary_vals_fn=expression_summary_vals_fn,
            cell_counts_fn=cell_counts_generator_fn,
        )

        cell_type_orderings = build_cell_orderings(cell_counts_cube_dir, cell_ordering_generator_fn)
        primary_filter_dimensions = build_precomputed_primary_filters()

        with tiledb.open(expression_summary_cube_dir, ctx=create_ctx()) as expression_summary_cube, tiledb.open(
            cell_counts_cube_dir, ctx=create_ctx()
        ) as cell_counts_cube:
            yield WmgSnapshot(
                snapshot_identifier=snapshot_name,
                expression_summary_cube=expression_summary_cube,
                expression_summary_default_cube=None,
                expression_summary_fmg_cube=None,
                marker_genes_cube=None,
                cell_counts_cube=cell_counts_cube,
                cell_type_orderings=cell_type_orderings,
                primary_filter_dimensions=primary_filter_dimensions,
                dataset_to_gene_ids=None,
            )


def build_cell_orderings(cell_counts_cube_dir_, cell_ordering_generator_fn) -> DataFrame:
    cell_type_orderings = []
    with tiledb.open(cell_counts_cube_dir_, ctx=create_ctx()) as cell_counts_cube:
        tissue_ontology_term_ids = cell_counts_cube.df[:]["tissue_ontology_term_id"].unique()
        for tissue_ontology_term_id in tissue_ontology_term_ids:
            cell_type_ontology_term_ids = sorted(
                cell_counts_cube.df[tissue_ontology_term_id]["cell_type_ontology_term_id"].unique()
            )
            ordering = cell_ordering_generator_fn(cell_type_ontology_term_ids)
            cell_type_orderings.append(
                pd.DataFrame(
                    data={
                        "tissue_ontology_term_id": [tissue_ontology_term_id] * len(cell_type_ontology_term_ids),
                        "cell_type_ontology_term_id": cell_type_ontology_term_ids,
                        "depth": list(islice(cycle([0, 1, 2]), len(cell_type_ontology_term_ids))),
                        "order": ordering,
                    }
                )
            )
    return pd.concat(cell_type_orderings)


def create_cubes(
    data_dir,
    dim_size: int = 3,
    dim_ontology_term_ids_generator_fn: Callable[[str, int], List[str]] = simple_ontology_terms_generator,
    exclude_logical_coord_fn: Callable[[List[str], Tuple], bool] = None,
    expression_summary_vals_fn: Callable[[List[Tuple]], Dict[str, List]] = random_expression_summary_values,
    cell_counts_fn: Callable[[List[Tuple]], List[int]] = random_cell_counts_values,
) -> Tuple[str, str]:
    coords, dim_values = build_coords(
        expression_summary_logical_dims, dim_size, dim_ontology_term_ids_generator_fn, exclude_logical_coord_fn
    )
    expression_summary_cube_dir = create_expression_summary_cube(
        data_dir, coords, dim_values, expression_summary_vals_fn=expression_summary_vals_fn
    )

    coords, dim_values = build_coords(
        cell_counts_logical_dims, dim_size, dim_ontology_term_ids_generator_fn, exclude_logical_coord_fn
    )
    cell_counts_cube_dir = create_cell_counts_cube(data_dir, coords, dim_values, cell_counts_fn=cell_counts_fn)

    return expression_summary_cube_dir, cell_counts_cube_dir


def create_cell_counts_cube(data_dir, coords, dim_values, cell_counts_fn: Callable[[List[Tuple]], List[int]]) -> str:
    cube_dir = f"{data_dir}/cell_counts"
    tiledb.Array.create(cube_dir, cell_counts_schema, overwrite=True)

    with tiledb.open(cube_dir, mode="w") as cube:
        logical_attr_values: Dict[str, list] = {
            "n_cells": cell_counts_fn(coords),
        }
        assert all([len(logical_attr_values[attr.name]) == len(coords) for attr in cell_counts_logical_attrs])

        physical_dim_values = dim_values[: len(cell_counts_indexed_dims)]
        physical_attr_values = {
            cell_counts_logical_dims[i]: dim_values[i]
            for i in range(len(cell_counts_indexed_dims), len(cell_counts_logical_dims))
        }

        physical_attr_values.update(logical_attr_values)
        cube[tuple(physical_dim_values)] = physical_attr_values

        return cube_dir


def create_expression_summary_cube(
    data_dir,
    coords,
    dim_values,
    expression_summary_vals_fn: Callable[[List[tuple]], Dict[str, List]] = random_expression_summary_values,
) -> str:
    cube_dir = f"{data_dir}/expression_summary"
    tiledb.Array.create(cube_dir, expression_summary_schema, overwrite=True)

    with tiledb.open(cube_dir, mode="w") as cube:
        logical_attr_values = expression_summary_vals_fn(coords)
        assert all([len(logical_attr_values[attr.name]) == len(coords) for attr in expression_summary_logical_attrs])

        physical_dim_values = dim_values[: len(expression_summary_indexed_dims)]
        physical_attr_values = {
            expression_summary_logical_dims[i]: dim_values[i]
            for i in range(len(expression_summary_indexed_dims), len(expression_summary_logical_dims))
        }
        physical_attr_values.update(logical_attr_values)
        cube[tuple(physical_dim_values)] = physical_attr_values

    return cube_dir


def build_coords(
    logical_dims,
    dim_size,
    dim_ontology_term_ids_generator_fn: Callable[[str, int], List[str]] = simple_ontology_terms_generator,
    exclude_coord_fn: Callable[[Tuple], bool] = None,
) -> Tuple[List[Tuple], List[List]]:
    n_dims = len(logical_dims)
    n_coords = dim_size**n_dims

    def dim_domain_values(i_dim: int, dim_size_: int) -> List[str]:
        dim_name = logical_dims[i_dim]
        domain_values = dim_ontology_term_ids_generator_fn(dim_name, dim_size_)
        assert len(set(domain_values)) == dim_size
        return domain_values

    all_dims_domain_values = [dim_domain_values(i_dim, dim_size) for i_dim in range(n_dims)]
    # create all possible coordinate values (dim_size ^ n_dims)
    dim_values = [
        [all_dims_domain_values[i_dim][(i_row // dim_size**i_dim) % dim_size] for i_row in range(n_coords)]
        for i_dim in range(n_dims)
    ]
    coords = list(zip(*dim_values))
    if exclude_coord_fn:
        Coord = namedtuple("Coord", logical_dims)
        coords: List[Tuple] = list(filterfalse(exclude_coord_fn, (Coord(*c) for c in coords)))
        dim_values: List[List] = [[coord_tuple[i_dim] for coord_tuple in coords] for i_dim in range(n_dims)]
    assert all([len(dim_values[i_dim]) == len(coords) for i_dim in range(n_dims)])
    return coords, dim_values


# CLI invocation for use by setup_dev_data.sh, to create a snapshot for Docker-based dev & test envs
if __name__ == "__main__":
    output_cube_dir = sys.argv[1]
    if not os.path.isdir(output_cube_dir):
        sys.exit(f"invalid dir {output_cube_dir} for cube")
    _, cell_counts_cube_dir = create_cubes(
        output_cube_dir,
        dim_size=4,
        dim_ontology_term_ids_generator_fn=semi_real_dimension_values_generator,
        exclude_logical_coord_fn=exclude_random_coords_75pct,
        expression_summary_vals_fn=random_expression_summary_values,
        cell_counts_fn=random_cell_counts_values,
    )
    cell_counts_df = build_cell_orderings(cell_counts_cube_dir, cell_ordering_generator_fn=forward_cell_type_ordering)
    cell_counts_df.to_json(os.path.join(output_cube_dir, CELL_TYPE_ORDERINGS_FILENAME), orient="records")
