import logging

import dask.distributed as dd

logger = logging.getLogger(__name__)

DASK_CLIENT = None


def start_dask_cluster():
    global DASK_CLIENT
    if DASK_CLIENT:
        return DASK_CLIENT
    logger.info("Starting Dask cluster")
    cluster = dd.LocalCluster()
    client = dd.Client(cluster)
    return client
