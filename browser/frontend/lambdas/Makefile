SHELL=/bin/bash -o pipefail
tests:=$(wildcard tests/test*.py)

ifndef APP_NAME
$(error Please set APP_NAME parameter using `make target APP_NAME=foo)
endif

.PHONY: deploy
deploy:
	$(MAKE) deploy -C chalice/$(APP_NAME)

.PHONY: destroy
destroy:
	$(MAKE) destroy -C chalice/$(APP_NAME)

.PHONY: clean
clean:
	$(MAKE) clean -C chalice/$(APP_NAME)
