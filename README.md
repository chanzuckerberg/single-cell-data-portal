# DCP Prototype
![Push Tests](https://github.com/HumanCellAtlas/dcp-prototype/workflows/Push%20Tests/badge.svg)
[![codecov](https://codecov.io/gh/chanzuckerberg/dcp-prototype/branch/master/graph/badge.svg)](https://codecov.io/gh/chanzuckerberg/dcp-prototype)

The Data Coordination Platform (DCP) stores, organizes, serves, and enables the exploration of cell atlases.

## Developer guide

### Configuration

Environment variables used for configuration:
* `DEPLOYMENT_STAGE` - set this value to target a specific deployment stage for deploying code and infrastructure.
    * `dev`, `staging`
* `AWS_PROFILE` - the AWS profile used to access and manage infrastructure
    * `single-cell-dev`, `single-cell-prod`

Configuring `awscli` credentials and profiles in `~/.aws`:
1. `pip install awscli`
1. Run `aws configure` or manually configure `czi-id` credentials:

    ```shell
    # ~/.aws/credentials example

    [czi-id]
    aws_access_key_id = ACCESS_KEY_ID
    aws_secret_access_key = SECRET_ACCESS_KEY
    ```

1.  [Configure](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-files.html)
    the `single-cell-dev` and/or `single-cell-prod` AWS profiles.

    ```shell
    # ~/.aws/config example

    [profile single-cell-dev]
    role_arn = arn:aws:iam::ACCOUNT_ID:role/poweruser
    source_profile = czi-id
    region = us-east-1

    [profile single-cell-prod]
    role_arn = arn:aws:iam::ACCOUNT_ID:role/poweruser
    source_profile = czi-id
    region = us-east-1
    ```

    Please contact #help-infra if you require access to `single-cell-dev` or `single-cell-prod`.

### Testing
Install dependencies `pip install -r requirements-dev.txt`

Unit tests:
* Run `make unit-test`

Functional tests:
* Set `DEPLOYMENT_STAGE` to a valid deployed environment (`dev`, `staging`)
* Run `make functional-test`

### Code formatting and linting

Commands and their uses:
* `make lint` - run code linter
* `make fmt` - run code auto-formatters

### Deploy
1. Set `DEPLOYMENT_STAGE`, and `AWS_PROFILE` in environment
1. Deploy Infra `make deploy -C infra`
1. [Deploy Backend](browser/backend/README.md#Deploy)
1. [Deploy Cloudfront-invalidator](browser/lambdas/README.md#Deploy)
1. [Deploy Frontend](browser/frontend/README.md#Deployment)
