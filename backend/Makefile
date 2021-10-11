.PHONY: db/migrate
db/migrate:
	PYTHONPATH=.. alembic -x db=${DEPLOYMENT_STAGE} -c=./config/database.ini upgrade head

.PHONY: db/remote_migration_init
db/remote_migration_init:
	# Run utility installation before invoking these commands.
	pip install awscli
	apt-get update && apt-get install -y postgresql-client jq

.PHONY: db/init_remote_dev
db/init_remote_dev: db/remote_migration_init
	$(eval DB_URI = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString' | jq -r .remote_dev_uri))
	# The create db script exits with a 1 status code if the db already exists.
	-@ echo "Creating DB..." && \
		cd .. && \
		python -m scripts.populate_db --create-schema --skip-populate --skip-recreate-db && \
		if [ ! -z "${DATA_LOAD_PATH}" ]; then \
			echo "Importing db snapshot from s3..." && \
			aws s3 cp ${DATA_LOAD_PATH} /tmp/db_snapshot.sql && \
			psql ${DB_URI}${REMOTE_DEV_PREFIX} < /tmp/db_snapshot.sql; \
		else \
			echo "Importing blank db snapshot..." && \
			psql ${DB_URI}${REMOTE_DEV_PREFIX} < scripts/db_snapshot.sql && \
			echo "Writing test data..." && \
			python -m scripts.populate_db --populate-data --skip-recreate-db; \
		fi

	PYTHONPATH=.. alembic -x db=${DEPLOYMENT_STAGE} -c=./config/database.ini upgrade head

.PHONY: db/delete_remote_dev
db/delete_remote_dev: db/remote_migration_init
	# Delete database.
	-@ echo "Deleting DB..." && \
		cd .. && \
		python -m scripts.populate_db --drop-db

db/rollback:
	PYTHONPATH=.. alembic -x db=${DEPLOYMENT_STAGE}  -c=./config/database.ini downgrade -1

db/new_migration:
	# Usage: make db/new_migration MESSAGE="purpose_of_migration"
	PYTHONPATH=.. alembic -c=./config/database.ini revision --message $(MESSAGE)

db/new_migration_auto:
	# Usage: make db/new_migration MESSAGE="purpose_of_migration"
	PYTHONPATH=.. alembic -c=./config/database.ini revision --autogenerate --message $(MESSAGE)

db/connect:
	# Assuming you've created a tunnel to the DB. Check the docs for information on how to do that
	$(eval DB_PW = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | split(":") | .[-1] | split("@") | .[0]'))
	PGPASSWORD=${DB_PW} psql --dbname corpora_${DEPLOYMENT_STAGE} --username corpora_${DEPLOYMENT_STAGE} --host 0.0.0.0 $(ARGS)

db/download:
    # Download the database to corpora_dev-<date>.sqlc
	$(eval DB_PW = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | split(":") | .[-1] | split("@") | .[0]'))
	$(eval OUTFILE = $(shell date +corpora_${DEPLOYMENT_STAGE}-%Y%m%d%H%M.sqlc))
	PGPASSWORD=${DB_PW} pg_dump -Fc --dbname=corpora_${DEPLOYMENT_STAGE} --file=database/${OUTFILE} --host 0.0.0.0 --username corpora_${DEPLOYMENT_STAGE}

db/import:
    # Imports corpora_dev.sqlc into the corpora_test database
	# Usage: make db/import FROM=corpora_dev-202102221309
	docker-compose exec database pg_restore --clean --no-owner --username corpora --dbname corpora /import/$(FROM).sqlc

db/import/schema:
    # Imports the corpora_dev.sqlc schema (schema ONLY) into the corpora_test database
	# Usage: DEPLOYMENT_STAGE=test make db/import/schema
	pg_restore --schema-only --clean --no-owner --dbname corpora_test corpora_$(DEPLOYMENT_STAGE).sqlc
	# Also import alembic schema version
	pg_restore --data-only --table=alembic_version --no-owner --dbname corpora_test corpora_$(DEPLOYMENT_STAGE).sqlc

db/dump_schema:
ifeq ($(DEPLOYMENT_STAGE),test)
	docker-compose exec database pg_dump --schema-only --dbname=corpora --username corpora
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

SSH_SERVER_ALIVE_INTERVAL_IN_SECONDS?=60
SSH_SERVER_ALIVE_COUNT_MAX?=60
db/tunnel:
	$(eval endpoint=$(shell aws rds describe-db-cluster-endpoints --db-cluster-identifier corpora-${DEPLOYMENT_STAGE}-corpora-api | jq '.DBClusterEndpoints[] | select(.EndpointType | contains("WRITER")) | .Endpoint'))
	ssh -N -o ServerAliveInterval=${SSH_SERVER_ALIVE_INTERVAL_IN_SECONDS} -o ServerAliveCountMax=${SSH_SERVER_ALIVE_COUNT_MAX} \
		-L 5432:${endpoint}:5432 bastion.${DEPLOYMENT_STAGE}.single-cell.czi.technology