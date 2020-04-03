SHELL:=/bin/bash

.PHONY: fmt
fmt:
	black .
	terraform fmt -recursive infra

.PHONY: lint
lint:
	flake8 dcp_prototype tests

.PHONY: unit-test
unit-test:
	PYTHONWARNINGS=ignore:ResourceWarning coverage run \
		--source=dcp_prototype/backend, \
		--omit=.coverage,venv,dcp_prototype/backend/browser/scripts,dcp_prototype/backend/browser/api \
		-m unittest discover \
		--start-directory tests/ \
		--top-level-directory . \
		--verbose
