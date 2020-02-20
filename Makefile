SHELL:=/bin/bash

.PHONY: fmt
fmt:
	black .

.PHONY: lint
lint:
	flake8 dcp_prototype tests

.PHONY: unit-test
unit-test:
	python -m pytest -s tests
