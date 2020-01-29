from uuid import uuid4

from .constants import ID_GENERATOR_PREFIX


def hca_accession_generator(entity):
    return ID_GENERATOR_PREFIX + "-" + entity + "-" + str(uuid4())
