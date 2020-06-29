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
		-m unittest discover --start-directory tests/unit/backend --top-level-directory . --verbose

.PHONY: functional-test
functional-test:
	python3 -m unittest discover --start-directory tests/functional --top-level-directory . --verbose
