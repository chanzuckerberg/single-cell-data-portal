import contextlib
import os
import sys
import tempfile
from itertools import filterfalse
from typing import List, Callable, Tuple

import numpy as np
import tiledb
from numpy.random import random, randint

from backend.corpora.common.corpora_orm import DbDataset, CollectionVisibility, DbCollection
from backend.corpora.common.entities import Collection
from backend.corpora.common.utils.db_session import db_session_manager
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
def create_temp_cube(
    dim_size=3,
    attr_vals_fn: Callable[[List[Tuple]], List] = random_attr_values,
    exclude_logical_coord_fn: Callable[[Tuple], bool] = None,
) -> None:
    with tempfile.TemporaryDirectory() as cube_dir:
        create_cube(cube_dir, dim_size, attr_vals_fn=attr_vals_fn, exclude_logical_coord_fn=exclude_logical_coord_fn)
        with tiledb.open(cube_dir) as cube:
            yield cube


def simple_ontology_terms_generator(dimension_name: str, n_terms: int) -> List[str]:
    return [f"{dimension_name}_{i}" for i in range(n_terms)]


def create_dataset(dataset_id_ordinal: int) -> str:
    coll_id = f"dataset_id_{dataset_id_ordinal}_coll_id"
    with db_session_manager() as session:

        if coll := Collection.get(session, (coll_id, CollectionVisibility.PUBLIC)):
            Collection.delete(coll)

        collection = DbCollection(
            id=coll_id,
            visibility=CollectionVisibility.PUBLIC.name,
            name=f"dataset_id_{dataset_id_ordinal}_coll_name",
            owner="owner",
        )
        session.add(collection)
        dataset = DbDataset(
            id=f"dataset_id_{dataset_id_ordinal}",
            name=f"dataset_name_{dataset_id_ordinal}",
            collection_id=coll_id,
            collection_visibility=CollectionVisibility.PUBLIC,
        )
        session.add(dataset)
        return dataset.id


def semi_real_dimension_values_generator(dimension_name: str, dim_size: int) -> List[str]:
    """
    Returns a set of ontology term ids, sampled from real ontologies. While these ontology terms
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
        return list(sorted(ontology_labels.gene_term_id_labels.keys()))[:dim_size]
    if dimension_name == "tissue_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("UBERON")][:dim_size]
    if dimension_name == "organism_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("NCBITaxon")][:dim_size]
    if dimension_name == "cell_type_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("CL")][:dim_size]
    if dimension_name == "dataset_id":
        return [create_dataset(i) for i in range(dim_size)]
    if dimension_name == "assay_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("EFO")][:dim_size]
    if dimension_name == "development_stage_ontology_term_id":
        return [
            term_id for term_id in deterministic_term_ids if term_id.startswith("Hsap") or term_id.startswith("MmusDev")
        ][:dim_size]
    if dimension_name == "disease_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("MONDO")][:dim_size]
    if dimension_name == "ethnicity_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("HANCESTRO")][:dim_size]
    if dimension_name == "sex_ontology_term_id":
        return [term_id for term_id in deterministic_term_ids if term_id.startswith("PATO")][:dim_size]
    raise AssertionError(f"unknown dimension name {dimension_name}")


def create_cube(
    cube_dir,
    dim_size=3,
    attr_vals_fn: Callable[[List[tuple]], List] = random_attr_values,
    dim_ontology_term_ids_generator_fn: Callable[[str, int], List[str]] = simple_ontology_terms_generator,
    exclude_logical_coord_fn: Callable[[Tuple], bool] = None,
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

        if exclude_logical_coord_fn:
            coord_tuples = list(filterfalse(exclude_logical_coord_fn, coord_tuples))
            coords = [[coord_tuple[i_dim] for coord_tuple in coord_tuples] for i_dim in range(n_dims)]

        assert all([len(coords[i_dim]) == len(coord_tuples) for i_dim in range(n_dims)])

        logical_attr_values = attr_vals_fn(coord_tuples)
        assert all([len(logical_attr_values[attr.name]) == len(coord_tuples) for attr in cube_logical_attrs])

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
    create_cube(output_cube_dir, dim_ontology_term_ids_generator_fn=semi_real_dimension_values_generator)
