from backend.layers.processing.process_cxg import ProcessCxg
from tests.unit.backend.fixtures.environment_setup import fixture_file_path

test_file = fixture_file_path("local_seurat_fail_oom.h5ad")  # replace this with the h5ad to profile

if __name__ == "__main__":
    """ """
    ProcessCxg(None, None, None).make_cxg(test_file)
