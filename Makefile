SHELL:=/bin/bash

black:
	black .

lint:
	flake8 --config=./.github/actions/flake8/.flake8 .
