# Cxg Schema Migrations

A Cell-by-Gene migration is a procedure that occurs when a new cellxgene schema is available. The schemas are defined [here](https://github.com/chanzuckerberg/single-cell-curation/tree/main/schema). In brief, a schema expresses requirements and restrictions on the structure, format, value types and interpretation of all single-cell datasets the Cell-by-Gene corpus. The migration procedure creates a new dataset, compliant with the latest schema release, for every dataset in our corpus.

This project is aimed at CZI Data Service engineers who are tasked with performing and monitoring the migration process. It assumes you have access to and have read the internal document outlining the migration procedure. This project is meant to assist those engineers with additional context and tools to monitor and diagnose in real-time, issues arising from migrations.

> 💡 **Note:** The migration process does issue a report with failures, but this is only available at the end of migration. Migration will take days.

## Pre-Requisites

- You have access to two AWS profiles, one for development environment (i.e. `dev`) and one for production environment (i.e. `prod`).
- You have install the AWS CLI

## Migration Environments

We perform a dry-run migration in the `dev` environment. The `dev` environment is a mirror of all the datasets available in `production` plus any additional datasets have been manually ingested into dev. For the remainder of this document, we assume `dev` is a 1:1 mirror of the data in `production`. The purpose of the `dev` migration is to catch, diagnose and squash all errors that arise during migration. Errors may arise from infrastructure issues or unanticipated migration logic issues.

Ultimately, the entire migration process must be applied to the `production` environment.

## Migration Procedure

The migration procedure consists of two stages (a) migration, and (b) ingestion.

The **migration stage** simply converts every dataset from the previous schema-version (i.e. 5.2.0) to the newest (i.e. 5.3.0) via a migration script. The migration script is tailor made by our team of curators and is available from the `cellxgene_schema` [Python package](https://pypi.org/project/cellxgene-schema/); pinned to the schema version to which it applies. The result is a dataset in our ecosystem, somewhere in storage.

The **ingestion stage** is responsible for bringing that migrated dataset into our registration ecosystem, persisted by a relational datase. It is the relational database that tracks and manages every dataset and their memership to both public and private collections. The same ingestion stage procedure handles both migrated datasets (ones that were previously in our corpus) and new datasets registered by our curators.

Both of these stages are encoded as state machines in AWS stepfunction service. These services and their codebase locations are listed in the table.

| stage     | step-function                        | python file                                   | entry point           |
| --------- | ------------------------------------ | --------------------------------------------- | --------------------- |
| migration | dp-dev-devstack-schema-migration-sfn | backend/layers/processing/schema_migration.py | SchemaMigrate.migrate |
| ingestion | dp-dev-devstack-sfn                  | backend/layers/processing/process.py          | ProcessMain.process   |

A simplified overview of these procedures are diagramed below:

<!-- TODO: Add Diagram -->

## Monitoring an ongoing migration

We kick off the migration by executing the migration stag step function. Each step in these step functions are executed as jobs, which are added to the same job-queue `schema_migration-STAGE_NAME` (where `STAGE_NAME` is one of `dev` or `prod`). The monitoring boils down to tracking all of the jobs enqueued by either the main migration stage step-function, or the downstream ingestion stage step-functions.

> 💡 **Note:** We currently don't tag a migration execution. Therefore we will use start times and end times to cull the relelvant jobs
