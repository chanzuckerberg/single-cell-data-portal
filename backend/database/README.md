# Data Portal Database Procedures

- General information about Alembic migrations can be found [here](https://alembic.sqlalchemy.org/en/latest/index.html).
- For database [make recipes](../Makefile)

## How to perform a database migration

1. From the top level directory (`corpora-data-portal`), `cd backend`.
1. Run `make db/new_migration MESSAGE="purpose_of_migration"` where `purpose_of_migration` is a short phrase describing why the database migration is occurring.
1. Head over to the file that was just created. The path will be something like `backend/database/versions/xxxxxxxxxxxx_purpose_of_migration.py`. Edit the `upgrade()` and `downgrade()` functions such that `upgrade()` contains the commands to perform the migration you would like and `downgrade()` contains the commands to undo it.
1. Rename the file by prepending the migration count to the filename (`xxx_purpose_of_migration` -> `03_xxx_purpose_of_migration`)
1. Update the `Revision ID` and the `revision` (used by alembic) to include the migration count
For example `Revision ID: a8cd0dc08805` becomes `Revision ID: 18_a8cd0dc08805` and `revision = "a8cd0dc08805"` becomes `revision = "18_a8cd0dc08805"` 
1. [Test your migration](#test-a-migration)
1. Check that [corpora.orm](../corpora/common/corpora_orm.py) matches up with your changes.
1. Once you've completed the changes, create a PR to get the functions reviewed. 
1. Once the PR is merged, you can run the migration.
1. [Connect to Remote RDS](#connect-to-remote-rds)
1. In a new terminal, complete the migration by running:
```shell
cd corpora-data-portal/backend
export CORPORA_LOCAL_DEV=1
export DEPLOYMENT_STAGE=test
make db/migrate
```

## How to autogenerate migration script

1. From the top level directory (`corpora-data-portal`), `cd backend`.
1. Make changes to [corpora_orm.py](../corpora/common/corpora_orm.py)
1. [Connect to Remote RDS](#connect-to-remote-rds)
1. Run Auto-migration:
```shell
cd corpora-data-portal/backend
export CORPORA_LOCAL_DEV=1
export DEPLOYMENT_STAGE=test
make db/new_migration_auto MESSAGE="purpose_of_migration"
```
See [What does Autogenerate Detect (and what does it not detect?)](https://alembic.sqlalchemy.org/en/latest/autogenerate.html#what-does-autogenerate-detect-and-what-does-it-not-detect).
5. Follow [How to perform a database migration](#how-to-perform-a-database-migration) starting from **step 3**.

### Test a Migration
The following steps will test that a migration script works on a local database using data downloaded from a deployed database. 

1. [Connect to Remote RDS](#connect-to-remote-rds)
2. Open a new terminal and download the remote database schema:
```shell
cd corpora-data-portal/backend
export DEPLOYMENT_STAGE=test
make db/download
```
This wll download the database into a file ending in *.sqlc*.
3. Close the tunnel to the remote database
4. Start the local database environment: 
```shell
cd corpora-data-portal
make local-start
export DEPLOYMENT_STAGE=test
```
5. Import the remote database schema into your local database:  
```shell
cd backend/
make db/import FROM=${}
```
where from is the name of the *.sqlc* file downloaded. For example 
```shell script
make db/import FROM=corpora_dev-202102221309
```
- Note: The file is stored under backend/database/file_name but the make command will look in import/file_name this is due to the way [the local paths are mapped to the docker container](https://github.com/chanzuckerberg/corpora-data-portal/blob/ffca067b9e4aea237fa2bd7c7a9cbc5813ebd449/docker-compose.yml#L13)

You may need to run this a few times, until there are no significant errors.
 - Note: `pg_restore: error: could not execute query: ERROR:  role "rdsadmin" does not exist` is not a significant error
6. Run the migration test:
```shell
make db/test_migration
``` 
If there are no differences then the test passed. If the test didn't pass, make adjustments to your migration script and restart from step 5. Repeat until there are no errors.

## Connect to Remote RDS
Enable local connection to the private RDS instance:

- Note: Since the default PostgreSQL port is `5432`, the above command will conflict with a local PostgreSQL instance.
To stop it run `make local-stop` from `./corpora-data-portal`


```shell
cd ./corpora-data-poral/backend
export DEPLOYMENT_STAGE=test
make db/tunnel
```

This command opens an SSH tunnel from `localhost:5432` to the RDS connection endpoint via the *bastion* server.
The local port `5432` is fixed and encoded in the DB connection string stored in the AWS Secret at
`corpora/backend/<DEPLOYMENT_STAGE>/database_local`.

