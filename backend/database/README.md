# Data Portal Database Procedures

- General information about Alembic migrations can be found [here](https://alembic.sqlalchemy.org/en/latest/index.html).
- For database [make recipes](../Makefile)
- See [Environment variables](../../README.md#environment-variables) for 
usage explanation of `DEPLOYMENT_STAGE`, `AWS_PROFILE` and `CORPORA_LOCAL_DEV`
- `$REPO_ROOT` - root directory where the `single-cell-data-portal` project is cloned (e.g. `~/PyCharmProjects/single-cell-data-portal`)

## How to perform a database migration

1. `cd $REPO_ROOT/backend`
1. Run `make db/new_migration MESSAGE="purpose_of_migration"` where `purpose_of_migration` is a short phrase describing why the database migration is occurring.
   This generates a file in `$REPO_ROOT/backend/database/versions/xxxxxxxxxxxx_purpose_of_migration.py`
1. In the generated file, edit the `upgrade()` and `downgrade()` functions such that `upgrade()` contains the Alembic DDL commands to perform the migration you would like and `downgrade()` contains the commands to undo it.
1. Rename the generated file by prepending the migration count to the filename (`xxx_purpose_of_migration.py` -> `03_xxx_purpose_of_migration.py`)
1. In the generated file, update the `Revision ID` and the `revision` (used by Alembic) to include the migration count.
For example `Revision ID: a8cd0dc08805` becomes `Revision ID: 18_a8cd0dc08805` and `revision = "a8cd0dc08805"` becomes `revision = "18_a8cd0dc08805"` 
1. [Test your migration](#test-a-migration)
1. Check that [corpora.orm.py](../corpora/common/corpora_orm.py) matches up with your changes.
1. Once you've completed the changes, create a PR to get the functions reviewed.
1. Once the PR is merged, you can run the migration.
1. [Connect to Remote RDS](#connect-to-remote-rds)
1. In a new terminal, complete the migration by running:
```shell
cd $REPO_ROOT/backend
CORPORA_LOCAL_DEV=1 DEPLOYMENT_STAGE=test make db/migrate
```

## How to autogenerate migration script

1. Make changes to the ORM class(es) in [corpora_orm.py](../corpora/common/corpora_orm.py)
2. [Connect to Remote RDS](#connect-to-remote-rds). Note, generally, you should be connecting to prod
   (`AWS_PROFILE=single-cell-prod DEPLOYMENT_STAGE=prod`) since we want to generate 
a migration from the database schema currently deployed in prod.
3. Autogenerate the migration using the steps below. `AWS_PROFILE` and `DEPLOYMENT_STAGE` should be the same values
used in the previous [Connect to Remote RDS](#connect-to-remote-rds) step. For details about Alembic's migration autogeneration, 
see [What does Autogenerate Detect (and what does it not detect?)](https://alembic.sqlalchemy.org/en/latest/autogenerate.html#what-does-autogenerate-detect-and-what-does-it-not-detect)
```shell
cd $REPO_ROOT/backend
AWS_PROFILE=single-cell-{dev,prod} DEPLOYMENT_STAGE={dev,staging,prod} CORPORA_LOCAL_DEV=1 make db/new_migration_auto MESSAGE="purpose_of_migration"
```
4. Follow [How to perform a database migration](#how-to-perform-a-database-migration) starting from **step 3** 
(i.e. editing the `upgrade()` and `downgrade()` functions).

### Test a Migration
The following steps will test that a migration script works on a local database using data downloaded from a deployed database. 

1. [Connect to Remote RDS](#connect-to-remote-rds)
2. Open a new terminal and using the same values for `AWS_PROFILE` and `DEPLOYMENT_STAGE`, download the remote database schema:
```shell
cd $REPO_ROOT/backend
AWS_PROFILE=single-cell-{dev,prod} DEPLOYMENT_STAGE={dev,staging,prod} make db/download
```
This will download the database into the `$REPO_ROOT/backend/database` directory in a file named `corpora_${DEPLOYMENT_STAGE}-<YYYYmmddHHMM>.sqlc`.

3. Close the tunnel to the remote database
4. Start the local database environment:
```shell
cd $REPO_ROOT
make local-start
```
5. Import the remote database schema into your local database:
```shell
cd $REPO_ROOT/backend
DEPLOYMENT_STAGE=test make db/import FROM=corpora_{dev,staging,prod}-<YYYYmmddHHMM>  # exclude the .sqlc extension
```
where the `FROM` parameter is the base name of the `.sqlc` file downloaded from the `make db/download` step above. For example 
```shell
make db/import FROM=corpora_prod-202102221309
```
- Note: The file is stored locally under `$REPO_ROOT/backend/database/corpora_${DEPLOYMENT_STAGE}-<YYYYmmddHHMM>.sqlc` 
but the `make db/import` command retrieves it from `/import/$(FROM).sqlc` due to the way [the local paths are mapped to the Docker container](https://github.com/chanzuckerberg/corpora-data-portal/blob/ffca067b9e4aea237fa2bd7c7a9cbc5813ebd449/docker-compose.yml#L13)

You may need to run this a few times, until there are no significant errors.
 - Note: `pg_restore: error: could not execute query: ERROR:  role "rdsadmin" does not exist` is not a significant error
6. Run the migration test:
```shell
AWS_PROFILE=single-cell-{dev,prod} DEPLOYMENT_STAGE=test make db/test_migration
``` 
This test will:
1. Dump the current schema (before)
1. Apply the migration (upgrade)
1. Rollback the migration (downgrade)
1. Dump the schema (after)
1. Compare the before vs after schemas. These should be identical if the database migration's `upgrade()` and `downgrade()` functions were implemented correctly.

If there are no differences then the test passed. If the test didn't pass, make adjustments to your migration script and restart from step 5. Repeat until there are no errors.

## Connect to Remote RDS
Enable local connection to the private RDS instance:

- Note: Since the default PostgreSQL port is `5432`, the above command will conflict with a local PostgreSQL instance.
To stop it run `make local-stop` from the `$REPO_ROOT` directory.


```shell
cd $REPO_ROOT/backend
AWS_PROFILE=single-cell-{dev,prod} DEPLOYMENT_STAGE={dev,staging,prod} make db/tunnel
```

This command opens an SSH tunnel from `localhost:5432` to the RDS connection endpoint via the *bastion* server.
The local port `5432` is fixed and encoded in the DB connection string stored in 
[AWS Secrets Manager](https://us-west-2.console.aws.amazon.com/secretsmanager/home?region=us-west-2#!/listSecrets/)
in the secret named `corpora/backend/${DEPLOYMENT_STAGE}/database_local`.

