"""Use to profile memory usage of process_cxg.py for a given h5ad file.
To run in the docker container use: $ make local-seurat-mem-profiler
"""
from backend.layers.processing.process_seurat import ProcessSeurat
from tests.unit.backend.fixtures.environment_setup import fixture_file_path

test_file = fixture_file_path("local_seurat_fail_oom.h5ad")  # replace this with the h5ad to profile

if __name__ == "__main__":
    ProcessSeurat(None, None, None).make_seurat(test_file)
