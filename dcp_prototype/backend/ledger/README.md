Requirements
```
export DEPLOYMENT_STAGE={environment_in_use}
export AWS_PROFILE={aws_profile_in_user}
```

## Database Management
```
# Delete and replace all tables
from sqlalchemy import create_engine
from dcp_prototype.backend.ledger.code.config.ledger_config import LedgerDbConfig

engine = create_engine(LedgerDbConfig().database_uri)
conn = engine.connect()
conn.execute('DROP SCHEMA public CASCADE;')
conn.execute('CREATE SCHEMA public;')
conn.execute('GRANT ALL ON SCHEMA public TO ledger_dev;')
conn.execute('GRANT ALL ON SCHEMA public TO public;')

# Run make db/migrate to re-populate all tables
```

## ORM Documentation
https://docs.sqlalchemy.org/en/13/orm/tutorial.html

Simple code samples are provided below. Please look at official documentation above for wide range of supported features and functions.

```
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
```