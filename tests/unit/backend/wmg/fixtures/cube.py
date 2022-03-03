import contextlib
import os
import sys
import tempfile
from typing import List, Callable, Dict

import numpy as np
import tiledb
from numpy.random import random, randint

from backend.wmg.data.schema import cube_logical_dims, schema, cube_indexed_dims, cube_logical_attrs


def random_attr_values(coords):
    return {
        "nnz": randint(size=len(coords), low=0, high=100),
        "n_cells": randint(size=len(coords), low=0, high=1000),
        "sum": random(size=len(coords)) * 10,
    }


def all_ones_attr_values(coords):
    attr_vals = {"nnz": np.ones(len(coords)), "n_cells": np.ones(len(coords)), "sum": np.ones(len(coords))}
    return attr_vals


@contextlib.contextmanager
def create_temp_cube(dim_size=3, attr_vals_fn: Callable[[List], Dict[str, list]] = random_attr_values) -> None:
    with tempfile.TemporaryDirectory() as cube_dir:
        create_cube(cube_dir, dim_size, attr_vals_fn)
        with tiledb.open(cube_dir) as cube:
            yield cube


def simple_ontology_terms_generator(dimension_name: str, n_terms: int) -> List[str]:
    return [f"{dimension_name}_{i}" for i in range(n_terms)]


# Functions that can be used to generate a set of valid ontology term ids, sampled from real ontologies. While these
# ontology terms from the real ontologies, they are not always ones that would be admissable w.r.t. to the cellxgene
# dataset schema. This is still useful to ensure that id-to-label mapping lookups will return a real label. Note that
# this implementation is wildly inefficient and in many cases does not return values that would be allowed in our
# real datasets, but it is good enough for test code.
def semi_real_ontology_terms_generator(dimension_name: str, n_terms: int) -> List[str]:
    # must import lazily
    import backend.wmg.data.ontology_labels as ontology_labels

    if ontology_labels.ontology_term_id_labels is None:
        ontology_labels.__load_ontologies()
    if ontology_labels.gene_term_id_labels is None:
        ontology_labels.__load_genes()

    deterministic_term_ids = sorted(ontology_labels.ontology_term_id_labels.keys())

    if dimension_name == "gene_ontology_term_id":
        return list(sorted(ontology_labels.gene_term_id_labels.keys()))[:n_terms]
    if dimension_name == "tissue_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("UBERON")][:n_terms]
    if dimension_name == "organism_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("NCBITaxon")][:n_terms]
    if dimension_name == "cell_type_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("CL")][:n_terms]
    if dimension_name == "dataset_id":
        return [f"dataset_id_{i}" for i in range(n_terms)]
    if dimension_name == "assay_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("EFO")][:n_terms]
    if dimension_name == "development_stage_ontology_term_id":
        return [
            term_id for term_id in deterministic_term_ids if term_id.startswith("Hsap") or term_id.startswith("MmusDev")
        ][:n_terms]
    if dimension_name == "disease_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("MONDO")][:n_terms]
    if dimension_name == "ethnicity_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("HANCESTRO")][:n_terms]
    if dimension_name == "sex_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("PATO")][:n_terms]
    raise AssertionError(f"unknown dimension name {dimension_name}")


def create_cube(
    cube_dir,
    dim_size=3,
    attr_vals_fn: Callable[[List[tuple]], List] = random_attr_values,
    dim_ontology_term_ids_generator_fn: Callable[[str, int], List[str]] = simple_ontology_terms_generator,
) -> None:
    if tiledb.array_exists(cube_dir):
        raise FileExistsError(cube_dir)

    tiledb.Array.create(cube_dir, schema)

    with tiledb.open(cube_dir, mode="w") as cube:
        n_dims = len(cube_logical_dims)
        n_coords = dim_size**n_dims

        def dim_domain_values(i_dim: int, dim_size_: int) -> List[str]:
            dim_name = cube_logical_dims[i_dim]
            domain_values = dim_ontology_term_ids_generator_fn(dim_name, dim_size_)
            assert len(set(domain_values)) == dim_size
            return domain_values

        all_dims_domain_values = [dim_domain_values(i_dim, dim_size) for i_dim in range(n_dims)]
        # create all possible coordinate values (dim_size ^ n_dims)
        coords = [
            [all_dims_domain_values[i_dim][(i_row // dim_size**i_dim) % dim_size] for i_row in range(n_coords)]
            for i_dim in range(n_dims)
        ]

        coord_tuples = list(zip(*coords))
        logical_attr_values = attr_vals_fn(coord_tuples)
        assert all([len(coords[i_dim]) == n_coords for i_dim in range(n_dims)])
        assert all([len(logical_attr_values[attr.name]) == n_coords for attr in cube_logical_attrs])

        physical_dim_values = coords[: len(cube_indexed_dims)]
        physical_attr_values = {
            cube_logical_dims[i]: coords[i] for i in range(len(cube_indexed_dims), len(cube_logical_dims))
        }
        physical_attr_values.update(logical_attr_values)
        cube[tuple(physical_dim_values)] = physical_attr_values


# CLI invocation for use by setup_dev_data.sh, to create a cube for Docker-based dev & test envs
if __name__ == "__main__":
    output_cube_dir = sys.argv[1]
    if not os.path.isdir(output_cube_dir):
        sys.exit(f"invalid dir {output_cube_dir} for cube")
    create_cube(output_cube_dir, dim_ontology_term_ids_generator_fn=semi_real_ontology_terms_generator)
