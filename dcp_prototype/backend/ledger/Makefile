.PHONY: db/migrate
db/migrate:
	alembic -x db=${DEPLOYMENT_STAGE} -c=./config/database.ini upgrade head

.PHONY: db/rollback
db/rollback:
	alembic -x db=${DEPLOYMENT_STAGE}  -c=./config/database.ini downgrade -1

.PHONY: db/new_migration
db/new_migration:
	# Usage: make db/new_migration MESSAGE="purpose_of_migration"
	alembic -c=./config/database.ini revision --message $(MESSAGE)

.PHONY: db/connect
db/connect:
	$(eval DATABASE_URI = $(shell aws secretsmanager get-secret-value --secret-id dcp/backend/ledger/${DEPLOYMENT_STAGE}/database --region us-east-1 | jq -r '.SecretString | fromjson.database_uri'))
	psql --dbname $(DATABASE_URI)

.PHONY: db/download
db/download:
	# Usage: make db/download FROM=dev    - downloads DB to ledger_dev-<date>.sqlc
	$(eval DATABASE_URI = $(shell aws secretsmanager get-secret-value --secret-id dcp/backend/ledger/${FROM}/database --region us-east-1 | jq -r '.SecretString | fromjson.database_uri'))
	$(eval OUTFILE = $(shell date +ledger_${FROM}-%Y%m%d%H%M.sqlc))
	pg_dump -Fc --dbname=$(DATABASE_URI) --file=$(OUTFILE)

.PHONY: db/import
db/import:
	# Usage: make db/import FROM=dev    - imports ledger_dev.sqlc into ledger_local
	pg_restore --clean --no-owner --dbname ledger_local ledger_$(FROM).sqlc

.PHONY: db/import/schema
db/import/schema:
	# Usage: DEPLOYMENT_STAGE=dev make db/import/schema  - imports ledger_dev.sqlc into ledger_local
	pg_restore --schema-only --clean --no-owner --dbname ledger_local ledger_$(DEPLOYMENT_STAGE).sqlc
	# Also import alembic schema version
	pg_restore --data-only --table=alembic_version --no-owner --dbname ledger_local ledger_$(DEPLOYMENT_STAGE).sqlc

.PHONY: db/dump_schema
db/dump_schema:
	pg_dump --schema-only --dbname=ledger_local

.PHONY: db/test_migration
db/test_migration:
	$(MAKE) db/dump_schema > /tmp/before
	$(MAKE) db/migrate
	$(MAKE) db/rollback
	$(MAKE) db/dump_schema > /tmp/after
	diff /tmp/{before,after} # No news is good news.
