SHELL:=/bin/bash

.PHONY: fmt
fmt:
	black .

.PHONY: lint
lint:
	flake8 backend tests

.PHONY: unit-test
unit-test:
	-docker run -d -p 5432:5432 --name test_db -e POSTGRES_PASSWORD=test_pw postgres
	PYTHONWARNINGS=ignore:ResourceWarning python3 -m coverage run \
		-m unittest discover --start-directory tests/unit/backend --top-level-directory . --verbose; \
	test_result=$$?; \
	$(MAKE) clean_test_db; \
	exit $$test_result

.PHONY: unittest
unittest:
	PYTHONWARNINGS=ignore:ResourceWarning python3 -m coverage run \
		-m unittest discover --start-directory tests/unit/backend --top-level-directory . --verbose; \
	test_result=$$?; \
	exit $$test_result

clean_test_db:
	-docker stop test_db
	-docker rm test_db

.PHONY: functional-test
functional-test:
	python3 -m unittest discover --start-directory tests/functional --top-level-directory . --verbose

.PHONY: local-database
local-database: clean_test_db
	docker run -d -p 5432:5432 --name test_db -e POSTGRES_PASSWORD=test_pw postgres
	python3 ./scripts/populate_db.py


.PHONY: local-backend
local-backend:
	$(MAKE) local-server -C ./backend/chalice/api_server

.PHONY: smoke-test-prod-build
smoke-test-prod-build:
	$(MAKE) smoke-test-prod-build -C ./frontend

.PHONY: smoke-test-with-local-backend
smoke-test-with-local-backend:
	$(MAKE) smoke-test-with-local-backend -C ./frontend

.PHONY: smoke-test-with-local-backend-ci
smoke-test-with-local-backend-ci:
	$(MAKE) smoke-test-with-local-backend-ci -C ./frontend

# ------ LOCAL-DEV ------
# Local-dev related commands are below this line for now.
oauth/pkcs12/certificate.pfx:
	# All calls to the openssl cli happen in the oidc-server-mock container.
	@echo "Generating certificates for local dev"
	docker run -ti -v $(PWD)/oauth/pkcs12:/tmp/certs --workdir /tmp/certs --rm=true --entrypoint bash soluto/oidc-server-mock ./generate_cert.sh
	@echo
	@echo "Installing generated certs into the local keychain requires sudo access:"
	sudo security add-trusted-cert -d -p ssl -k /Library/Keychains/System.keychain oauth/pkcs12/server.crt
	docker run -ti -v $(PWD)/oauth/pkcs12:/tmp/certs --workdir /tmp/certs --rm=true --entrypoint bash soluto/oidc-server-mock ./generate_pfx.sh

.PHONY: dev-init
dev-init: oauth/pkcs12/certificate.pfx
	docker-compose up -d
	docker-compose exec backend pip3 install awscli
	docker-compose exec backend /corpora-data-portal/scripts/setup_dev_data.sh

.PHONY: dev-status
dev-status:
	docker ps -a | grep --color=no -e 'CONTAINER\|corpora-data-portal'

.PHONY: dev-sync
dev-sync: dev-init
	docker-compose up --build -d

.PHONY: dev-start
dev-start:
	docker-compose up -d

.PHONY: dev-stop
dev-stop:
	docker-compose stop

.PHONY: dev-clean
dev-clean:
	-if [ -f ./oauth/pkcs12/server.crt ] ; then \
	    export CERT=$$(docker run -ti -v $(PWD)/oauth/pkcs12:/tmp/certs --workdir /tmp/certs --rm=true --entrypoint "" soluto/oidc-server-mock bash -c "openssl x509 -in server.crt -outform DER | sha1sum | cut -d ' ' -f 1"); \
	    echo ""; \
	    echo "Removing this certificate requires sudo access"; \
	    sudo security delete-certificate -Z $${CERT} /Library/Keychains/System.keychain; \
	fi;
	-rm -rf ./oauth/pkcs12/server*
	-rm -rf ./oauth/pkcs12/certificate*
	docker-compose rm -sf
	-docker volume rm corpora-data-portal_database
	-docker volume rm corpora-data-portal_localstack

# make dev-logs or make dev-logs CONTAINER=backend
.PHONY: dev-logs
dev-logs:
	docker-compose logs -f $(CONTAINER)

.PHONY: dev-shell
dev-shell:
	docker-compose exec $(CONTAINER) bash

.PHONY: dev-unit-test
dev-unit-test:
	docker-compose exec backend bash -c "cd /corpora-data-portal && make unittest"

.PHONY: dev-functional-test
dev-functional-test:
	docker-compose exec backend bash -c "cd /corpora-data-portal && export DEPLOYMENT_STAGE=test && make functional-test"

.PHONY: dev-smoke-test
dev-smoke-test:
	docker-compose exec frontend make smoke-test-with-local-dev
