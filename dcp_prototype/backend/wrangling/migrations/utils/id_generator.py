from uuid import uuid4

from .constants import ID_GENERATOR_PREFIX


def hca_accession_generator(entity):
    """
    Generates a unique accession number for each DCP entity using the entity type as
    part of the accession.
    """

    return ID_GENERATOR_PREFIX + "-" + entity + "-" + str(uuid4())


def hca_accession_transformer(entity, uuid):
    """
    Generates an accession number retaining a given uuid using the standard pattern for
    all HCA accessions.
    """
    return ID_GENERATOR_PREFIX + "-" + entity + "-" + str(uuid)
