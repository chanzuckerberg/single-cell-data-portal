#!/usr/bin/env python

import os
import sys

import click

pkg_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))  # noqa
sys.path.insert(0, pkg_root)  # noqa

from backend.corpora.common.utils.db_session import db_session_manager
from backend.corpora.common.entities.dataset import Dataset


@click.command()
@click.argument('uuid')
def delete_dataset(uuid):
    """Delete a dataset from Cellxgene"""
    click.echo(uuid)
    with db_session_manager() as session:
        dataset = Dataset.get(session, uuid)
        print(dataset.to_dict())
        # dataset.dataset_and_asset_deletion()
        # dataset.delete()


if __name__ == '__main__':
    delete_dataset()
