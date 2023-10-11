include ../../common.mk

# Bundles the application so that it can be dockerized and run in ECS

cwd := $(shell pwd)

.PHONY: build-ecs
build-ecs:
	$(call copy_client_assets,../../client/build,../../server/common/web) ; \
	cp ../../client/build/index.html ../../server/common/web/templates/index.html ; \
	cp -r ../../client/build/* ../../server/common/web/. ;
