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

.PHONY: package
package:
	$(MAKE) package -C api
