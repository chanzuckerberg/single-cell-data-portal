# LOCAL DEVELOPMENT ENVIRONMENT WITHOUT DOCKER

## Setting up local python virtual environment for WMG development

### Install python packages needed for WMG

1. Make sure that you have successfully completed the [prerequisite installation steps](./README.md#pre-requisite-installations-and-setups)
1. Create a python virtual environment using your preferred method. Some options are [`venv`](https://realpython.com/python-virtual-environments-a-primer/) and [`miniconda`](https://conda.io/projects/conda/en/stable/user-guide/install/macos.html#install-macos-silent)

1. **Install pygraphviz** -
   There is currently a problem with installing `pygraphviz` on Mac OSX. This package is needed for the WMG pipeline. The installation problem and the solution is explained [here](https://github.com/pygraphviz/pygraphviz/issues/11#issuecomment-1380458670)

   - `brew install graphviz`
   - `pip install --global-option=build_ext --global-option="-I$(brew --prefix graphviz)/include/" --global-option="-L$(brew --prefix graphviz)/lib/" pygraphviz==1.11`

1. Install base packages - `pip install -r requirements.txt`
1. Install packages for WMG api - `pip install -r backend/api_server/requirements.txt`
1. Install packages for WMG pipeline - `pip install -r backend/wmg/pipeline/requirements.txt`

### Run unit tests

1. Run unit tests for WMG api - `pytest -v tests/unit/backend/wmg`
1. Run unit tests for WMG pipeline - `pytest -v tests/unit/wmg_processing`

### Run functional tests

Run functional tests for WMG api against the `dev` environment.
**NOTE**: `dev` environment is a remote environment. These functional tests run locally against a backend in a remote environment called `dev`.

1. `AWS_PROFILE=single-cell-dev DEPLOYMENT_STAGE=dev pytest -v tests/functional/backend/wmg/test_wmg_api.py`
