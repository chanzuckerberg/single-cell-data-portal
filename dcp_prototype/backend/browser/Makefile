SHELL=/bin/bash -o pipefail
tests:=$(wildcard tests/test*.py)

.PHONY: deploy
deploy:
	$(MAKE) deploy -C api

.PHONY: destroy
destroy:
	$(MAKE) destroy -C api

.PHONY: clean
clean:
	$(MAKE) clean -C api

.PHONY: test
test: package $(tests)
	coverage combine
	rm -f .coverage.*

.PHONY: package
package:
	$(MAKE) package -C api

$(tests): %.py : lint
	coverage run -p --source=browser $*.py -v

.PHONY: unit-test
unit-test: package
	python -m unittest discover -p test_*api.py -t . -v
