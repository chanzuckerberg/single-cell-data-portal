import time

from backend.wmg.pipeline.summary_cubes.snapshot_builder import WmgSnapshotBuilder

ORGANISM_INFO = [
    {"label": "homo_sapiens", "id": "NCBITaxon:9606"},
    {"label": "mus_musculus", "id": "NCBITaxon:10090"},
]


def run():
    snapshot_id = str(int(time.time()))
    for organism in ORGANISM_INFO:
        builder = WmgSnapshotBuilder(corpus_path=snapshot_id, organismInfo=organism)
        builder.run_pipeline()
