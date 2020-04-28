SHELL=/bin/bash -o pipefail
tests:=$(wildcard tests/test*.py)

ifndef app
$(error Please set APP_NAME parameter using `make target APP_NAME=foo)
endif

.PHONY: deploy
deploy:
	$(MAKE) deploy -C $APP_NAME

.PHONY: destroy
destroy:
	$(MAKE) destroy -C $APP_NAME

.PHONY: clean
clean:
	$(MAKE) clean -C $APP_NAME
