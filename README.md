# Corpora Data Portal
![Push Tests](https://github.com/chanzuckerberg/corpora-data-portal/workflows/Push%20Tests/badge.svg)
[![codecov](https://codecov.io/gh/chanzuckerberg/corpora-data-portal/branch/main/graph/badge.svg)](https://codecov.io/gh/chanzuckerberg/corpora-data-portal)
[![Maintainability](https://api.codeclimate.com/v1/badges/9416c313de4d0457a5cc/maintainability)](https://codeclimate.com/github/chanzuckerberg/corpora-data-portal/maintainability)

The Corpora Data Portal (CDP) enables the publication, discovery and exploration of interoperable
single-cell datasets. Data contributors can upload, review and publish datasets for private or
public use. Via the portal, data consumers are able to discover, download and connect data to visualization tools
such as [cellxgene](https://chanzuckerberg.github.io/cellxgene/posts/cellxgene_cziscience_com) to perform further
analysis. The goal of the CDP is to catalyze distributed collaboration of single-cell research by providing a large,
well-labeled repository of interoperable datasets.

## Developers

### Pre-requisites
1. Python >= 3.6
1. [install moreutils](https://joeyh.name/code/moreutils/)
1. [Install and configure awscli](docs/awscli.md)
1. [Configure ssh access](https://github.com/chanzuckerberg/single-cell-infra#ssh)
1. [Set development flags](docs/dev_flags.md)

### Environment variables
| Name | Description | Values |
|------|-------------|--------|
|`DEPLOYMENT_STAGE`|Specifies an app deployment stage for tasks such as deployments and functional tests.|`dev`, `staging`, `prod`|
|`AWS_PROFILE`|Specifies the profile used to interact with AWS resources via awscli.|`single-cell-dev`, `single-cell-prod`|
|`CORPORA_LOCAL_DEV`|If this variable is set to any value, the Corpora app will operate in local development mode.|Any|

### Commands
| Command | Description | Notes |
|---------|-------------|-------|
|`make fmt`|Auto-format codebase using [black](https://pypi.org/project/black/).|This should be run before merging in any changes.|
|`make lint`|Perform lint checks on codebase using [flake8](https://flake8.pycqa.org/en/latest/).|This should be run before merging in any changes.|
|`make unit-test`|Run all unit tests.||
|`make functional-tests`|Run all functional tests.|These tests run against a deployed environment which is selected by the value of `DEPLOYMENT_STAGE`.|

### Development
1. [Backend](backend/chalice/api_server/README.md#Development)
1. [Frontend](frontend/README.md#Development)

### Deployment
1. Set `DEPLOYMENT_STAGE` and `AWS_PROFILE` according to the environment to be deployed.
1. [Deploy Backend](backend/chalice/api_server/README.md#Deploy)
1. [Deploy Cloudfront-invalidator](backend/chalice/cloudfront_invalidator/README.md#Deploy)
1. [Deploy Frontend](frontend/README.md#Deployment)


### Running unittests
1. Set `AWS_PROFILE`
1. Build the chalice app for the api-server `$ make package -C ./backend/chalice/api_server`
1. Run the tests `$ make unit-test`
