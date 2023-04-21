import psutil
import tiledb

from backend.common.utils.math_utils import GB, MB


def virtual_memory_size(vm_fraction: float) -> int:
    mem_size_bytes = psutil.virtual_memory().total
    fractional_mem_bytes = int(vm_fraction * mem_size_bytes)
    return fractional_mem_bytes // MB * MB  # round down to MB boundary


def consolidation_buffer_size(vm_fraction: float) -> int:
    # consolidation buffer heuristic to prevent thrashing: total_mem/io_concurrency_level, rounded to GB
    io_concurrency_level = int(tiledb.Config()["sm.io_concurrency_level"])
    buffer_size = int(virtual_memory_size(vm_fraction) / io_concurrency_level) + (GB - 1)
    return buffer_size // GB * GB  # round down to GB boundary
