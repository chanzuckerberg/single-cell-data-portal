db/migrate:
	alembic -x db=${DEPLOYMENT_STAGE} -c=./config/database.ini upgrade head

db/rollback:
	alembic -x db=${DEPLOYMENT_STAGE}  -c=./config/database.ini downgrade -1

db/new_migration:
	# Usage: make db/new_migration MESSAGE="purpose_of_migration"
	alembic -c=./config/database.ini revision --message $(MESSAGE)

db/connect:
	$(eval DATABASE_URI = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | fromjson.database_uri'))
	psql --dbname=$(DATABASE_URI)

db/download:
    # Download the database to corpora_dev-<date>.sqlc
	# Usage: make db/download FROM=dev
	$(eval DATABASE_URI = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/${FROM}/database --region us-west-2 | jq -r '.SecretString | fromjson.database_uri'))
	$(eval OUTFILE = $(shell date +corpora_${FROM}-%Y%m%d%H%M.sqlc))
	pg_dump -Fc --dbname=$(DATABASE_URI) --file=$(OUTFILE)

db/import:
    # Imports corpora_dev.sqlc into the corpora_test database
	# Usage: make db/import FROM=dev
	pg_restore --clean --no-owner --dbname corpora_test corpora_$(FROM).sqlc

db/import/schema:
    # Imports the corpora_dev.sqlc schema (schema ONLY) into the corpora_test database
	# Usage: DEPLOYMENT_STAGE=dev make db/import/schema
	pg_restore --schema-only --clean --no-owner --dbname corpora_test corpora_$(DEPLOYMENT_STAGE).sqlc
	# Also import alembic schema version
	pg_restore --data-only --table=alembic_version --no-owner --dbname corpora_test corpora_$(DEPLOYMENT_STAGE).sqlc

db/dump_schema:
	pg_dump --schema-only --dbname=corpora_test

db/test_migration:
	$(MAKE) db/dump_schema > /tmp/before
	$(MAKE) db/migrate
	$(MAKE) db/rollback
	$(MAKE) db/dump_schema > /tmp/after
	diff /tmp/{before,after} # No news is good news.
