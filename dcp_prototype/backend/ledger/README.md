Ledger

ORM Documentation
https://docs.sqlalchemy.org/en/13/orm/tutorial.html

Simple code samples are provided below. Please look at official documentation above for wide range of supported features and functions.

{code}
from dcp_prototype.backend.ledger.code.common.ledger_orm import DBSessionMaker, BiosamplePrep, LibraryPrepProtocol

session_maker = DBSessionMaker()
session = session_maker.session()

# Query for record
record = session.query(BiosamplePrep).filter(BiosamplePrep.id == 'biosample_prep_id').one_or_none()
record = db.query(BiosamplePrep).get('biosample_prep_id')

# Insert Record
biosample_prep = BiosamplePrep(id='123', category='123')
session.add(biosample_prep)
session.commit()

# Delete Record
record = db.query(BiosamplePrep).get('biosample_prep_id')
session.delete(record)
session.commit()