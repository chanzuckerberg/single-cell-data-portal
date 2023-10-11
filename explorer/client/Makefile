include ../common.mk

CXG_CONFIG := $(if $(CXG_CONFIG),$(CXG_CONFIG),./__tests__/e2e/test_config.yaml)
SERVER_PORT_FILE := ../.test_server_port.txt

# Packaging
.PHONY: clean
clean:
	rm -rf node_modules
	rm -f __tests__/screenshots/*.png

.PHONY: ci
ci:
	npm ci

.PHONY: install
install:
	npm install

.PHONY: build
WEBPACK_CONFIG ?= configuration/webpack/webpack.config.prod.js
build:
	npm run build $(WEBPACK_CONFIG)

# start local frontend server
# nb: if using a new/updated dataset, must first start local backend/API server (`launch_dev_server.sh <dataset_path>`)
.PHONY: start-frontend
start-frontend:
	@-pbcopy < ../.test_base_url.txt;
	@if [ -e $(SERVER_PORT_FILE) ]; then \
		port=$$(cat $(SERVER_PORT_FILE)); \
	fi; \
	export CXG_SERVER_PORT=$${port:-$$CXG_SERVER_PORT}; \
	node server/development.js

# start an instance of the Explorer server and run the end-to-end tests
.PHONY: smoke-test
smoke-test:
	start_server_and_test \
		'../launch_dev_server.sh --config-file $(CXG_CONFIG) --port $(CXG_SERVER_PORT)' \
		$(CXG_SERVER_PORT) \
		'CXG_URL_BASE="http://localhost:$(CXG_SERVER_PORT)" npm run e2e -- --verbose --no-cache false'

.PHONY: unit-test
unit-test:
	node node_modules/jest/bin/jest.js --testPathIgnorePatterns e2e

# pass remaining commands through to npm run
%:
	npm run $(*)
