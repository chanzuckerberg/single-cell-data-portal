import logging
import time
from contextlib import contextmanager

try:
    from server_timing import Timing as ServerTiming
except ImportError:

    class ServerTiming:
        @staticmethod
        @contextmanager
        def time(key):  # noqa
            yield


logging.basicConfig(level=logging.INFO)


@contextmanager
def log_time_taken(description: str = "Code block"):
    start_time = time.time()
    try:
        with ServerTiming.time(description):
            yield
    finally:
        end_time = time.time()
        elapsed_time = end_time - start_time
        logging.info(dict(type="METRIC", details=dict(description=description, time=elapsed_time, unit="seconds")))
