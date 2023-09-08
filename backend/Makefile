.PHONY: db/migrate
db/migrate:
	PYTHONPATH=.. alembic -x db=${DEPLOYMENT_STAGE} -c=./database/database.ini upgrade head

.PHONY: db/remote_migration_init
db/remote_migration_init:
	# Run utility installation before invoking these commands.
	pip install awscli -r ../requirements-backend.txt
	apt-get update && apt-get install -y postgresql-client jq

.PHONY: db/init_remote_dev
db/init_remote_dev: db/remote_migration_init
	$(eval DB_URI = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString' | jq -r .remote_dev_uri))
	# The create db script exits with a 1 status code if the db already exists.
	-@ echo "Creating DB..." && \
		cd .. && \
		python3 -m scripts.populate_db --create-schema --skip-populate --skip-recreate-db && \
		if [ ! -z "${DATA_LOAD_PATH}" ]; then \
			echo "Importing db snapshot from s3..." && \
			aws s3 cp ${DATA_LOAD_PATH} /tmp/db_snapshot.sql && \
			psql ${DB_URI}${REMOTE_DEV_PREFIX} < /tmp/db_snapshot.sql; \
		else \
			echo "Importing blank db snapshot..." && \
			psql ${DB_URI}${REMOTE_DEV_PREFIX} < scripts/db_snapshot.sql && \
			echo "Writing test data..." && \
			python3 -m scripts.populate_db --populate-data --skip-recreate-db; \
		fi

	PYTHONPATH=.. alembic -x db=${DEPLOYMENT_STAGE} -c=./database/database.ini upgrade head

.PHONY: db/delete_remote_dev
db/delete_remote_dev: db/remote_migration_init
	# Delete database.
	-@ echo "Deleting DB..." && \
		cd .. && \
		python3 -m scripts.populate_db --drop-db

db/rollback:
	PYTHONPATH=.. alembic -x db=${DEPLOYMENT_STAGE}  -c=./database/database.ini downgrade -1

db/new_migration:
	# Usage: make db/new_migration MESSAGE="purpose_of_migration"
	PYTHONPATH=.. alembic -c=./database/database.ini revision --message "$(MESSAGE)"

db/new_migration_auto:
	# Usage: make db/new_migration_auto MESSAGE="purpose_of_migration"
	PYTHONPATH=.. alembic -c=./database/database.ini revision --autogenerate --message "$(MESSAGE)"

db/check:
	# Check if the database needs to be migrated due to changes in the schema.
	PYTHONPATH=.. alembic -c=./database/database.ini check

# interactive mode usage: AWS_PROFILE=single-cell-dev DEPLOYMENT_STAGE=dev make db/connect
# ARGS usage: AWS_PROFILE=single-cell-dev DEPLOYMENT_STAGE=dev make db/connect ARGS="-c \"select * from dataset_artifact where filetype='CXG'\""
db/connect:
ifeq ($(DEPLOYMENT_STAGE), rdev)
	$(MAKE) db/connect_internal DB_NAME=/${STACK} DB_USER=dataportal
else
	$(MAKE) db/connect_internal DB_NAME=corpora_${DEPLOYMENT_STAGE} DB_USER=corpora_${DEPLOYMENT_STAGE}
endif

db/connect_internal:
	$(eval DB_PW = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | match(":([^:]*)@").captures[0].string'))
	$(MAKE) db/tunnel/up
	PGOPTIONS='-csearch_path=persistence_schema' PGPASSWORD=${DB_PW} psql --dbname ${DB_NAME} --username ${DB_USER} --host 0.0.0.0 $(ARGS)
	$(MAKE) db/tunnel/down

db/console: db/connect # alias

db/dump:
  # Dump the DEPLOYMENT_STAGE database to OUTFILE
	$(eval DB_PW = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | match(":([^:]*)@").captures[0].string'))
	$(MAKE) db/tunnel/up
	PGPASSWORD=${DB_PW} pg_dump -Fc --dbname=corpora_${DEPLOYMENT_STAGE} --file=$(OUTFILE) --host 0.0.0.0 --username corpora_${DEPLOYMENT_STAGE}
	$(MAKE) db/tunnel/down

db/local/load-data:
	# Loads corpora_dev.sqlc into the local Docker env corpora database
	# Usage: make db/local/load-data INFILE=<file>
	$(eval DB_PW = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/test/database --region us-west-2 | jq -r '.SecretString | match(":([^:]*)@").captures[0].string'))
	PGPASSWORD=${DB_PW} pg_restore --clean --no-owner --host 0.0.0.0 --username corpora --dbname corpora $(INFILE)

db/local/load-schema:
    # Imports the corpora_dev.sqlc schema (schema ONLY) into the corpora_test database
	# Usage: make db/local/load-schema INFILE=<file>
    $(eval DB_PW = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/test/database --region us-west-2 | jq -r '.SecretString | match(":([^:]*)@").captures[0].string'))
	PGPASSWORD=${DB_PW} pg_restore --schema-only --schema=persistence-schema --clean --no-owner --host 0.0.0.0 --username corpora --dbname corpora $(INFILE)
	# Also import alembic schema version
	PGPASSWORD=${DB_PW} pg_restore --data-only --table=alembic_version --no-owner --host 0.0.0.0 --username corpora --dbname corpora $(INFILE)

db/dump_schema:
ifeq ($(DEPLOYMENT_STAGE),test)
	docker-compose exec database pg_dump --schema-only --dbname=corpora --username corpora
else
	$(eval DB_PW = $(shell aws secretsmanager get-secret-value --secret-id corpora/backend/${DEPLOYMENT_STAGE}/database --region us-west-2 | jq -r '.SecretString | match(":([^:]*)@").captures[0].string'))
	$(MAKE) db/tunnel/up
	PGPASSWORD=${DB_PW} pg_dump --schema-only --dbname corpora_${DEPLOYMENT_STAGE} --username corpora_${DEPLOYMENT_STAGE} --host 0.0.0.0
	$(MAKE) db/tunnel/down
endif

db/test_migration:
	$(MAKE) db/dump_schema > /tmp/before_migration
	$(MAKE) db/migrate
	$(MAKE) db/dump_schema > /tmp/after_migration
	$(MAKE) db/rollback
	$(MAKE) db/dump_schema > /tmp/after_rollback
	diff /tmp/{before_migration,after_rollback} # No news is good news.


SSH_SERVER_ALIVE_INTERVAL_IN_SECONDS?=60
SSH_SERVER_ALIVE_COUNT_MAX?=60
SSH_SOCKET=/tmp/data-portal-ssh-db-tunnel-socket-${DEPLOYMENT_STAGE}
ifeq ($(DEPLOYMENT_STAGE), rdev)
	SSH_BASTION_HOST=bastion.dev.single-cell.czi.technology
	CLUSTER_NAME=dataportal-rdev-happy

else
	SSH_BASTION_HOST=bastion.${DEPLOYMENT_STAGE}.single-cell.czi.technology
	CLUSTER_NAME=corpora-${DEPLOYMENT_STAGE}-corpora-api

endif
# TODO:
# - add db/tunnel as a dependency for all targets so that a tunnel is automatically opened if not already
db/tunnel/up:
	$(eval endpoint=$(shell aws rds describe-db-cluster-endpoints --db-cluster-identifier ${CLUSTER_NAME} | jq -r '.DBClusterEndpoints[] | select(.EndpointType | contains("WRITER")) | .Endpoint'))
	ssh -f -T -N -M -S $(SSH_SOCKET)\
		-o ServerAliveInterval=${SSH_SERVER_ALIVE_INTERVAL_IN_SECONDS} -o ServerAliveCountMax=${SSH_SERVER_ALIVE_COUNT_MAX} \
		-o ExitOnForwardFailure=yes \
		-L 5432:${endpoint}:5432 $(SSH_BASTION_HOST)

db/tunnel: db/tunnel/up # alias for backwards compatibility

db/tunnel/down:
	ssh -S $(SSH_SOCKET) -O exit $(SSH_BASTION_HOST) || true

SRC_ENV := prod
mirror_env_data:
	# Mirrors the SRC_ENV env's AWS RDS database and S3 data to
	# DEST_ENV. Defaults to prod->dev.
	#
	# If DEST_ENV = 'rdev':
	# - STACK arg must be included accordingly
	# - A comma-separated list of Collection id uuids may be passed as the COLLECTIONS arg to determine which Dataset
	# assets are mirrored from SRC_ENV s3 to the specified rdev STACK -- pass any non-null value for the DATA arg
	# (e.g., DATA=1) to copy the real datafiles, otherwise a sample Dataset is used to populate all relevant s3 asset
	# keys for all corresponding Datasets (faster). Covers all Dataset versions in each Collection. Ignored for mirrors
	# where DEST_ENV != rdev.
	#
	# THIS IS DESTRUCTIVE for the DEST_ENV env! The SRC_ENV env will
	# never be modified, but the DEST_ENV env's data will be replaced.
	#
	# Usage: make mirror_env_data [SRC_ENV={prod|staging|dev}] [DEST_ENV={dev|staging|rdev}] [STACK=<rdev_stack_name> [COLLECTIONS=<uuid1,uuid2,...> [DATA=1]]]
	scripts/mirror_env_data.sh $(SRC_ENV) $(DEST_ENV) $(STACK) $(COLLECTIONS) $(DATA)
