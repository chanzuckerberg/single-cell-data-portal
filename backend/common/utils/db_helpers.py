import logging

from sqlalchemy.orm import Session

from backend.common.corpora_orm import Base, DbDatasetProcessingStatus

logger = logging.getLogger(__name__)


def clone(model: Base, primary_key: dict = None, **kwargs) -> Base:
    """Clone an arbitrary sqlalchemy model object with new primary keys.
    https://stackoverflow.com/questions/28871406/how-to-clone-a-sqlalchemy-db-object-with-new-primary-key

    :param model: The SQLAlchemy model to clone
    :param primary_key: The new primary key values. If None a new primary key is generated
    :param kwargs: Updates the columns in the cloned model.
    :return: a clone of the original model with any kwargs passed in and the new primary key.
    """
    table = model.__table__
    non_pk_columns = [key for key in table.columns if key not in table.primary_key]
    data = {column: getattr(model, column) for column in non_pk_columns}
    data.update(kwargs)
    if primary_key:
        data.update(primary_key)
    return model.__class__(**data)


def processing_status_updater(session: Session, uuid: str, updates: dict):
    """
    Updates the DbDatasetProcessingStatus object with the specified statuses
    Note: The uuid parameter is the ID of the dataset_processing_status row, _NOT_ the dataset!
    :param session: active SQLAlchemy session
    :param uuid: ID of the dataset_processing_status to update
    :param updates: dictionary of statuses to update
    :return: None
    """
    session.query(DbDatasetProcessingStatus).filter(DbDatasetProcessingStatus.id == uuid).update(updates)
    session.commit()
    logger.debug("updating status", updates)
