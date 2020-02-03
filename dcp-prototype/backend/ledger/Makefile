db/migrate:
	alembic -x db=${DEPLOYMENT_STAGE} -c=./config/database.ini upgrade head

db/rollback:
	alembic -x db=${DEPLOYMENT_STAGE}  -c=./config/database.ini downgrade -1

db/new_migration:
	# Usage: make db/new_migration MESSAGE="purpose_of_migration"
	alembic -c=./config/database.ini revision --message $(MESSAGE)

db/connect:
	$(eval DATABASE_URI = $(shell aws secretsmanager get-secret-value --secret-id dcp/upload/${DEPLOYMENT_STAGE}/database --region us-east-1 | jq -r '.SecretString | fromjson.database_uri'))
	psql --dbname $(DATABASE_URI)

db/download:
	# Usage: make db/download FROM=dev    - downloads DB to upload_dev-<date>.sqlc
	$(eval DATABASE_URI = $(shell aws secretsmanager get-secret-value --secret-id dcp/upload/${FROM}/database --region us-east-1 | jq -r '.SecretString | fromjson.database_uri'))
	$(eval OUTFILE = $(shell date +upload_${FROM}-%Y%m%d%H%M.sqlc))
	pg_dump -Fc --dbname=$(DATABASE_URI) --file=$(OUTFILE)

db/import:
	# Usage: make db/import FROM=dev    - imports upload_dev.sqlc into upload_local
	pg_restore --clean --no-owner --dbname upload_local upload_$(FROM).sqlc

db/import/schema:
	# Usage: DEPLOYMENT_STAGE=dev make db/import/schema  - imports upload_dev.sqlc into upload_local
	pg_restore --schema-only --clean --no-owner --dbname upload_local upload_$(DEPLOYMENT_STAGE).sqlc
	# Also import alembic schema version
	pg_restore --data-only --table=alembic_version --no-owner --dbname upload_local upload_$(DEPLOYMENT_STAGE).sqlc

db/dump_schema:
	pg_dump --schema-only --dbname=upload_local

db/test_migration:
	$(MAKE) db/dump_schema > /tmp/before
	$(MAKE) db/migrate
	$(MAKE) db/rollback
	$(MAKE) db/dump_schema > /tmp/after
	diff /tmp/{before,after} # No news is good news.