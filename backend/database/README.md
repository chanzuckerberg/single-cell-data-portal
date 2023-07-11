# Data Portal Database Procedures

- General information about Alembic migrations can be found [here](https://alembic.sqlalchemy.org/en/latest/index.html).
- For database [make recipes](../Makefile)
- See [Environment variables](../../README.md#environment-variables) for
  usage explanation of `DEPLOYMENT_STAGE`, `AWS_PROFILE` and `CORPORA_LOCAL_DEV`
- `$REPO_ROOT` - root directory where the `single-cell-data-portal` project is cloned (e.g. `~/PyCharmProjects/single-cell-data-portal`)

## How to perform a database migration

1. `cd $REPO_ROOT/backend`
2. Run `make db/new_migration MESSAGE="purpose_of_migration"` where `purpose_of_migration` is a short phrase describing why the database migration is occurring.
   This generates a file in `$REPO_ROOT/backend/database/versions/xxxxxxxxxxxx_purpose_of_migration.py`
3. In the generated file, edit the `upgrade()` and `downgrade()` functions such that `upgrade()` contains the Alembic DDL commands to perform the migration you would like and `downgrade()` contains the commands to undo it.
4. Rename the generated file by prepending the migration count to the filename (`xxx_purpose_of_migration.py` -> `03_xxx_purpose_of_migration.py`)
5. In the generated file, update the `Revision ID` and the `revision` (used by Alembic) to include the migration count.
   For example `Revision ID: a8cd0dc08805` becomes `Revision ID: 18_a8cd0dc08805` and `revision = "a8cd0dc08805"` becomes `revision = "18_a8cd0dc08805"`
6. [Test your migration](#test-a-migration)
7. Check that [orm.py](../layers/persistence/orm.py) matches up with your changes.
8. Once you've completed the changes, create a PR to get the functions reviewed.
9. Rdev compatibility:
   1. If your migration requires changes to the seed data file for rdev, make those changes in a **new** seed data file with your migration hash in the name.
      ```
      s3://env-rdev-dataportal/database/seed_data_##_a1b2c3d4e5f6.sql
      ```
      where `##` is the migration number and `a1b2c3d4e5f6` is the migration hash, both taken from the migration file under `./versions/`. Do not delete older seed files from s3; the most recent version(s) will still be in active use by other developers before they incorporate your migration into their feature branches.
   2. In your PR, change the `data_load_path` local variable [in the ecs-stack module](../../.happy/terraform/modules/ecs-stack/main.tf) to reflect the new dump file above that you have written to s3.
10. Once the PR is merged, migrations will be run as part of the deployment process to each env.

While migrations should run automatically as part of our deployment process, if a manual migration is required:

1. [Connect to Remote RDS](#connect-to-remote-rds) to single-cell-dev
2. In a new terminal, complete the migration in the single-cell-dev test env by running:

```shell
cd $REPO_ROOT/backend
DEPLOYMENT_STAGE=test make db/migrate
```

## How to autogenerate migration script

1. Make changes to the ORM class(es) in [orm.py](../layers/persistence/orm.py)
2. [Connect to Remote RDS](#connect-to-remote-rds). Note, generally, you would be connecting to prod
   (`AWS_PROFILE=single-cell-prod DEPLOYMENT_STAGE=prod`) since we want to generate
   a migration from the database schema currently deployed in prod. However, if there are migrations haven't been
   deployed to prod yet, you would connect to staging here.
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

1. Open a new terminal and using the same values for `AWS_PROFILE` and `DEPLOYMENT_STAGE`, download the remote dev database schema:

```shell
cd $REPO_ROOT/backend
AWS_PROFILE=single-cell-{dev,prod} DEPLOYMENT_STAGE={dev,staging,prod} make db/dump OUTFILE=corpora_dev.sqlc
```

This will download the database to `$REPO_ROOT/backend/corpora_dev.sqlc`.

2. The tunnel to dev should close automatically (but worth verifying `ps ax | grep ssh`)
3. Start the local database environment:

```shell
cd $REPO_ROOT
make local-start
```

4. Import the remote database schema into your local database:

```shell
cd $REPO_ROOT/backend
make db/local/load-schema INFILE=corpora_dev.sqlc
```

where the `INFILE` parameter is the base name of the `.sqlc` file downloaded from the `make db/dump` step above. For example

```shell
make db/local/load-schema INFILE=corpora_dev.sqlc
```

You may need to run this a few times, until there are no significant errors.

- Note: `pg_restore: error: could not execute query: ERROR: role "rdsadmin" does not exist` is not a significant error

5. Run the migration test:

```shell
AWS_PROFILE=single-cell-{dev,prod} DEPLOYMENT_STAGE=test make db/test_migration
```

This test will:

1. Dump the current schema (before)
1. Apply the migration (upgrade)
1. Rollback the migration (downgrade)
1. Dump the schema (after)
1. Compare the before vs after schemas. These should be identical if the database migration's `upgrade()` and `downgrade()` functions were implemented correctly.

If there are no differences then the test passed. If the test didn't pass, make adjustments to your migration script and restart from step 4. Repeat until there are no errors.

## Connect to Remote RDS

Enable local connection to the private RDS instance:

- Note: Since the default PostgreSQL port is `5432`, the above command will conflict with a local PostgreSQL instance.
  To stop it run `make local-stop` from the `$REPO_ROOT` directory.

```shell
cd $REPO_ROOT/backend
AWS_PROFILE=single-cell-{dev,prod} DEPLOYMENT_STAGE={dev,staging,prod} make db/tunnel
```

This command opens an SSH tunnel from `localhost:5432` to the RDS connection endpoint via the _bastion_ server.
The local port `5432` is fixed and encoded in the DB connection string stored in
[AWS Secrets Manager](https://us-west-2.console.aws.amazon.com/secretsmanager/home?region=us-west-2#!/listSecrets/)
in the secret named `corpora/backend/${DEPLOYMENT_STAGE}/database`.
