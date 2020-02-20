SHELL=/bin/bash -o pipefail
tests:=$(wildcard tests/test*.py)

deploy:
	$(MAKE) deploy -C api

destroy:
	$(MAKE) destroy -C api

clean:
	$(MAKE) clean -C api

test: package $(tests)
	coverage combine
	rm -f .coverage.*

package:
	$(MAKE) package -C api

$(tests): %.py : lint
	coverage run -p --source=browser $*.py -v

unittest: package
	python -m unittest discover -p test_*api.py -t . -v

.PHONY: deploy package init destroy clean lint test
