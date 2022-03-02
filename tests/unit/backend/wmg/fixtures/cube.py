import contextlib
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


def create_cube(cube_dir, dim_size=3, attr_vals_fn: Callable[[List[tuple]], List] = random_attr_values) -> None:
    if tiledb.array_exists(cube_dir):
        raise FileExistsError(cube_dir)

    tiledb.Array.create(cube_dir, schema)

    with tiledb.open(cube_dir, mode="w") as cube:
        n_dims = len(cube_logical_dims)
        n_coords = dim_size**n_dims

        def dim_domain_values(i_dim: int, dim_size_: int) -> List[str]:
            domain_values = [f"{cube_logical_dims[i_dim]}_{i}" for i in range(dim_size_)]
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
