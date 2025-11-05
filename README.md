# CZ CELLxGENE Discover

![Push Tests](https://github.com/chanzuckerberg/single-cell-data-portal/workflows/Push%20Tests/badge.svg)
[![codecov](https://codecov.io/gh/chanzuckerberg/single-cell-data-portal/branch/main/graph/badge.svg?token=iIXh8Rw0CH)](https://codecov.io/gh/chanzuckerberg/single-cell-data-portal)

CELLxGENE Discover enables the publication, discovery and exploration of interoperable single-cell datasets. Data contributors can upload, review and publish datasets for private or public use. Via Discover, data consumers are able to discover, download and connect data to visualization tools such as [CELLxGENE Explorer](https://github.com/chanzuckerberg/cellxgene-documentation/blob/main/README.md) to perform further analysis. The goal of Discover is to catalyze distributed collaboration of single-cell research by providing a large, well-labeled repository of interoperable datasets.

## Code of Conduct

This project adheres to the Contributor Covenant [code of conduct](https://github.com/chanzuckerberg/.github/blob/master/CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to [opensource@chanzuckerberg.com](mailto:opensource@chanzuckerberg.com).

## Developer Guidelines

### FE Setup

Follow instructions in [/frontend README](./frontend/README.md)

### Pre-requisite installations and setups

**Note**: Before you begin to install any Python packages, make sure you have activated your Python virtual environment.

1. Install pre-commit: `pre-commit install` or check doc [here](https://pre-commit.com/)
1. Set up your machine to be able to work with AWS using the instructions [here](https://czi.atlassian.net/wiki/spaces/DC/pages/332892073/Getting+started+with+AWS). Please ensure to follow the step 3 `AWS CLI access` instructions all the way to the bottom so that you are also set up for SSH access. When you run the final command that requires the team's infra repo, use `single-cell-infra`.
1. [install jq](https://stedolan.github.io/jq/download/). If brew is installed run `brew install jq`.
1. [install libpq](https://www.timescale.com/blog/how-to-install-psql-on-mac-ubuntu-debian-windows/). If brew is installed run ` brew install libpq`. We need this
   tool to invoke the `psql` commandline utility.
1. Install chamber. For running functional tests below, you will need to install Chamber on your machine. Chamber is a tool for reading secrets stored in AWS Secret Store and Parameter Store. On Linux, go to https://github.com/segmentio/chamber/releases to download the latest version >= 2.9.0, and add it somewhere on your path. On Mac, run `brew install chamber`.
1. Install `uv` for Python dependency management (see below).

   ```bash
   # Using pip
   pip install uv

   # Using pipx (recommended)
   pipx install uv

   # Using Homebrew (macOS)
   brew install uv
   ```

### Python Dependency Management

This project uses `uv` for dependency management to ensure reproducible builds across all Docker images.

#### Dependency File Structure

- **`requirements.in`** - Source files containing direct dependencies (may include version constraints)
- **`requirements.txt`** - Auto-generated locked files with exact versions and hashes for all dependencies (direct + transitive)

All dependency files are located in the `/python_dependencies/` directory, organized by service:

#### Adding or Updating Dependencies

1. **Edit the source file** - Modify the appropriate `.in` file to add, update, or remove dependencies
2. **Regenerate the locked file** - Use the make update-dependencies recipe to generate the locked .txt file:
3. **Commit both files** - Always commit both the `.in` and `.txt` files together:

### Development Quickstart

**Note:** Make sure you are running your Python virtual environment before going through the development guides.

Once you have run the pre-requisite sets, you are ready to begin developing for CELLxGENE Discover. As you start to change code, you may want to deploy a test instance of Discover so that you can check to see how your changes perform. We have two ways to deploy your changes:

1. **Creating a local python-only development environment.** This local development environment is a minimal environment that is suitable for running the automated test suite quickly. It also allows you to run specific test cases which is not possible in the dockerized local deployment environment (see below) - the `make` commands invoke docker commands to run test suites in their entirety. Moreover, creating an environment without docker also makes for much faster test runs. To set up a python-only local development environment, see [this guide](DEV_ENV_WITHOUT_DOCKER.md)

1. **Creating a local deployment environment.** This environment will be entirely hosted on your own machine. It relies upon Docker to run both CELLxGENE Discover servers, Discover unit tests, and infrastructure service dependencies (AWS, Postgres, OIDC). The environment will be initialized with a small amount of dummy data. This environment is great to have up and running while you are actively developing. See [this guide](DEV_ENV.md) for instructions on how to set up a local deployment.

1. **Creating a remote deployment.** This environment creates a lightweight replica of CELLxGENE Discover, hosted by AWS, and provides a more realistic test bed to test your changes before either sending them to a PR or try them out with a cross-functional partner. It takes a longer time to deploy your changes to a remote development environment which is why the local deployment is preferred until your changes are ready for broader review. See [this guide](https://docs.google.com/document/d/1nynGcBS_TA55qlQo9WjINGkcMnE_xIBz7-inmop2bqo/edit#) for instructions on how to set up an rDev environment.

### Common Commands

| Command                                                                                 | Description                                                                          | Notes                                             |
| --------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------ | ------------------------------------------------- |
| `make fmt`                                                                              | Auto-format codebase using [black](https://pypi.org/project/black/).                 | This should be run before merging in any changes. |
| `make lint`                                                                             | Perform lint checks on codebase using [flake8](https://flake8.pycqa.org/en/latest/). | This should be run before merging in any changes. |
| `make unit-test`                                                                        | Run all unit tests.                                                                  |                                                   |
| `DEPLOYMENT_STAGE=<deployed_env> python3 -m unittest discover tests/functional/backend` | Run all functional tests.                                                            |                                                   |

### Environment variables

Environment variables are set using the command `export <name>=<value>`. For example, `export DEPLOYMENT_STAGE=dev`. These environment variables typically need to be set before you are able to set up your environments (i.e. local, rDev) and before you are able to successfully run any test suite.

| Name                | Description                                                                                                                                                                                   | Values                                |
| ------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------- |
| `DEPLOYMENT_STAGE`  | Specifies an app deployment stage for tasks such as deployments and functional tests. The `test` value implies local Docker development environment (and should probably be renamed `local`). | `test`, `dev`, `staging`, `prod`      |
| `AWS_PROFILE`       | Specifies the profile used to interact with AWS resources via the awscli.                                                                                                                     | `single-cell-dev`, `single-cell-prod` |
| `CORPORA_LOCAL_DEV` | Flag: If this variable is set to any value, the app will look for the database on **localhost:5432** and will use the aws secret `corpora/backend/\${DEPLOYMENT_STAGE}/database`.             | Any                                   |

### Database Procedures

If you need to make a change to the CELLxGENE Discover database, see [CELLxGENE Discover Database Procedures](backend/database/README.md).

### Running Unit Tests

1. Set `AWS_PROFILE`.
1. Ensure that you have set up your local development environment per the instructions above and run `make local-init` to launch a local dev environment.
1. Run the tests using the command `$ make unit-test`.

### Running the WMG pipeline and endpoints locally in an interactive environment

Please take a look at the tutorial notebooks in `example_dev_notebooks` for examples on how to run the WMG pipeline and endpoints locally for development purposes.

### Troubleshooting Unit Tests

1. If your unit tests crash due to an `Error 137`, that means the docker containers currently running are using up more memory than what the docker application
   has allocated. First run `$ make local-stop` to kill all docker containers and rerun the tests. If that doesn't work either, you may need to increase the
   memory allocation by going into _settings_ pane of your desktop docker application. A reasonable allocation of memory is `16GB`.

### Running Functional Tests

1. Set `AWS_PROFILE`.
2. Set `DEPLOYMENT_STAGE` to deployed environment you want to run tests, as written locally, against (dev, staging, or prod)
3. Run a specific suite of tests using `DEPLOYMENT_STAGE=<deployed_env> python3 -m unittest <path_to_functional_test>`. For example, `DEPLOYMENT_STAGE=dev python3 -m unittest tests/functional/backend/corpora/test_revisions.py`
4. Run all functional tests by using `DEPLOYMENT_STAGE=<deployed_env> python3 -m unittest discover tests/functional/backend`

### Running e2e Tests

Follow instructions in [e2e tests README](frontend/tests/README.md)

## Trademarks

CZ CELLXGENE, CZ CELLXGENE DISCOVER, and CZ CELLXGENE ANNOTATE are trademarks of the Chan Zuckerberg Initiative. All rights reserved.

Use, reuse, modification, and re-distribution of the source code in this repository is subject to the terms of the applicable open source [license](LICENSE). However, that license does not grant permission to use the trademarks without separate, express permission from the Chan Zuckerberg Initiative.
