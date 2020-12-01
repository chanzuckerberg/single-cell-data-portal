.PHONY: db/create
db/create:
	# TODO(mbarrien): Figure out which directory init_dev_db should live in
	python3 ../scripts/init_dev_db.py

.PHONY: db/migrate
db/migrate:
	alembic -x db=${DEPLOYMENT_STAGE} -c=./config/database.ini upgrade head

.PHONY: db/init_remote_dev
db/init_remote_dev: db/create db/migrate

db/rollback:
	alembic -x db=${DEPLOYMENT_STAGE}  -c=./config/database.ini downgrade -1

db/new_migration:
	# Usage: make db/new_migration MESSAGE="purpose_of_migration"
	alembic -c=./config/database.ini revision --message $(MESSAGE)

db/new_migration_auto:
	# Usage: make db/new_migration MESSAGE="purpose_of_migration"
	alembic -c=./config/database.ini revision --autogenerate --message $(MESSAGE)

db/connect:
	# Assuming you've created a tunnel to the DB. Check the docs for information on how to do that
	$(eval DB_PW = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | split(":") | .[-1] | split("@") | .[0]'))
	PGPASSWORD=${DB_PW} psql --dbname corpora_${DEPLOYMENT_STAGE} --username corpora_${DEPLOYMENT_STAGE} --host 0.0.0.0

db/download:
    # Download the database to corpora_dev-<date>.sqlc
	$(eval DB_PW = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | split(":") | .[-1] | split("@") | .[0]'))
	$(eval OUTFILE = $(shell date +corpora_${DEPLOYMENT_STAGE}-%Y%m%d%H%M.sqlc))
	PGPASSWORD=${DB_PW} pg_dump -Fc --dbname=corpora_${DEPLOYMENT_STAGE} --file=${OUTFILE} --host 0.0.0.0 --username corpora_${DEPLOYMENT_STAGE}

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
ifeq ($(DEPLOYMENT_STAGE),"test")
	pg_dump --schema-only --dbname=corpora_test
else
	$(eval DB_PW = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | split(":") | .[-1] | split("@") | .[0]'))
	PGPASSWORD=${DB_PW} pg_dump --schema-only --dbname corpora_${DEPLOYMENT_STAGE} --username corpora_${DEPLOYMENT_STAGE} --host 0.0.0.0
endif


db/test_migration:
	$(MAKE) db/dump_schema > /tmp/before
	$(MAKE) db/migrate
	$(MAKE) db/rollback
	$(MAKE) db/dump_schema > /tmp/after
	diff /tmp/{before,after} # No news is good news.
