## Data Portal Database Procedures

General information about Alembic migrations can be found [here](https://alembic.sqlalchemy.org/en/latest/index.html).

### How to perform a database migration

1. From the top level directory (`corpora-data-portal`), `cd backend`.
2. Run `make db/new_migration MESSAGE="purpose_of_migration"` where `purpose_of_migration` is a short phrase describing why the database migration is occurring.
3. Head over to the file that was just created. The path will be something like `backend/database/versions/xxxxxxxxxxxx_purpose_of_migration.py`. Edit the `upgrade()` and `downgrade()` functions such that `upgrade()` contains the commands to perform the migration you would like and `downgrade()` contains the commands to undo it.
4. Once you've completed the changes (and at this point, you can create a PR to get the functions reviewed), test your migration. Ensure that you are [connected](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/backend/chalice/api_server/README.md#development) to the database you'd like to migrate (usually `TEST`). Run `DEPLOYMENT_STAGE=test make db/migrate`.
5. If no errors are outputted, complete the migration by running `DEPLOYMENT_STAGE=dev make db/migrate`. Ensure that tunneling is set up and `CORPORA_LOCAL_DEV` is set to true.

### How to autogenerate migration script

1. From the top level directory (`corpora-data-portal`), `cd backend`.
1. Makes changes to ./corpora/common/corpora_orm.py
1. Set DEPLOYMENT_STAGE=dev, CORPORA_LOCAL_DEV=1
1. Run `make db/connect`
1. Follow [How to perform a database migration](#how_to_perform_a_database_migration) instructions except use `make db/new_migration_auto MESSAGE="purpose_of_migration"`. See [What does Autogenerate Detect (and what does it not detect?)](https://alembic.sqlalchemy.org/en/latest/autogenerate.html#what-does-autogenerate-detect-and-what-does-it-not-detect).
