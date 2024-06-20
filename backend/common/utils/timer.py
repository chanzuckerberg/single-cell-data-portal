import logging
import time
from contextlib import contextmanager
from typing import Any

import numpy as np

logging.basicConfig(level=logging.INFO)

Accumulators: dict[str, Any] = dict()


@contextmanager
def log_time_taken(
    description: str = "Code block",
):
    start_time = time.time()
    try:
        yield
    finally:
        end_time = time.time()
        elapsed_time = end_time - start_time
        logging.info(dict(type="METRIC", details=dict(description=description, time=elapsed_time, unit="seconds")))


class TimeAccumulator:
    def __init__(self, group: str):
        self.group = group
        self.records = []

    @contextmanager
    def time(
        self,
        description: str = "Code block",
    ):
        start_time = time.time()
        try:
            yield
        finally:
            end_time = time.time()
            logging.info(
                dict(
                    type="METRIC",
                    details=dict(description=description, time=end_time - start_time, unit="seconds", group=self.group),
                )
            )
            self.records.append(dict(description=description, elapsed=end_time - start_time))

    def todict(self):
        records = self.records
        records.sort(key=lambda x: x["elapsed"])
        times = [record["elapsed"] for record in records]
        self.total_time = sum(times)
        self.max_time = max(times)
        self.min_time = min(times)
        self.median_time = np.median(times)
        self._95th_percentile, self._99th_percentile = np.percentile(times, [95, 99])
        return dict(
            group=self.group,
            total_time=self.total_time,
            max_time=self.max_time,
            min_time=self.min_time,
            median_time=self.median_time,
            _95th_percentile=self._95th_percentile,
            _99th_percentile=self._99th_percentile,
            unit="seconds",
        )


@contextmanager
def log_accumulative_time_taken(group: str) -> TimeAccumulator:
    accumulator = TimeAccumulator(group)
    try:
        yield accumulator
    finally:
        logging.info(dict(type="METRIC", details=accumulator.todict()))
