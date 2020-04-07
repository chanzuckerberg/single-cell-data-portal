# DCP Prototype
![Push Tests](https://github.com/HumanCellAtlas/dcp-prototype/workflows/Push%20Tests/badge.svg)
![Code Coverage](https://codecov.io/gh/humancellatlas/dcp-prototype/branch/master/graph/badge.svg)

The Data Coordination Platform (DCP) stores, organizes, serves, and enables the exploration of cell atlases.

## Developer guide

### Configuration

Environment variables used for configuration:
* `DEPLOYMENT_STAGE` - set this value to target a specific deployment stage for deploying code and infrastructure.

### Testing
Set the `DEPLOYMENT_STAGE` environment variable to `test`

Install dependencies `pip install -r requirements-dev.txt`

Commands and their uses:
* `make unit-test` - runs all unit tests
* `make functional-test` - runs all functional tests

### Code formatting and linting

Commands and their uses:
* `make lint` - run code linter
* `make fmt` - run code auto-formatters
