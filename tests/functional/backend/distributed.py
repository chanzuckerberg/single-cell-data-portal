import json
from typing import Callable

from filelock import FileLock


def distributed_singleton(tmp_path_factory, worker_id: str, func: Callable) -> dict:
    """
    This function wraps a pytest fixture so it is only instantiated once and shared across all workers in a distributed
    test run.
    """
    if worker_id == "master":
        # not executing with multiple workers, just produce the data and let
        # pytest's fixture caching do its job
        return func()
    # get the temp directory shared by all workers
    root_tmp_dir = tmp_path_factory.getbasetemp().parent

    fn = root_tmp_dir.joinpath(func.__name__ + ".json")
    with FileLock(str(fn) + ".lock"):
        if fn.is_file():
            data = json.loads(fn.read_text())
        else:
            data = func()
            fn.write_text(json.dumps(data))
    return data
