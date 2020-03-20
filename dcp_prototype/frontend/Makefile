SHELL:=/bin/bash
export BROWSER_S3=dcp-static-site-test-699936264352

ifndef DEPLOYMENT_STAGE
$(error Please set the DEPLOYMENT_STAGE environment before deploying)
endif

deploy: init
	npm run build
	aws s3 sync ./public s3://$(BROWSER_S3)/

.PHONY: init
init:
	npm install

.PHONY: clean
clean:
	npm run clean
