# LOCAL DEVELOPMENT ENVIRONMENT WITHOUT DOCKER

## Setting up local python virtual environment for WMG development

### Install python packages needed for WMG

1. Make sure that you have successfully completed the [prerequisite installation steps](./README.md#pre-requisite-installations-and-setups)
1. _Create and activate_ a python virtual environment using your preferred method. Some options are [`venv`](https://realpython.com/python-virtual-environments-a-primer/) and [`miniconda`](https://conda.io/projects/conda/en/stable/user-guide/install/macos.html#install-macos-silent). Here is a [cheatsheet of commands](https://conda.io/projects/conda/en/stable/user-guide/cheatsheet.html) to manage your virtual environment with `miniconda`.

   **NOTE:** _As of this writing on 06/12/2023, you should create your virtual environment with_ `python 3.10` _or greater_

   - If you have installed `miniconda`, here is a list of basic virtual environment management commands:

     ```
     $ conda create -n <env-name> python=3.10 # creates virtual environment

     $ conda activate <env-name> # activates virtual environment

     $ conda deactivate # deactivate virtual environment
     ```

   - If you are using `venv`, here is a list of basic virtual environment management commands:

     ```
     $ cd single-cell-data-portal # root of this code repo

     $ python3 -m venv venv/ # create virtual environment

     $ source venv/bin/activate # activate virtual environment

     $ deactivate # deactivate virtual environment
     ```

1. **Install pygraphviz** -
   There is currently a problem with installing `pygraphviz` on Mac OSX. This package is needed for the WMG pipeline. The installation problem and the solution is explained [here](https://github.com/pygraphviz/pygraphviz/issues/11#issuecomment-1380458670)

   - `brew install graphviz`
   - `pip install --global-option=build_ext --global-option="-I$(brew --prefix graphviz)/include/" --global-option="-L$(brew --prefix graphviz)/lib/" pygraphviz==1.11`

1. Install packages for WMG api - `pip install -r requirements-backend.txt`
1. Install packages for WMG pipeline - `pip install -r requirements-wmg-pipeline.txt`

### Run unit tests for WMG

1. Run unit tests for WMG api - `pytest -v tests/unit/backend/wmg`
1. Run unit tests for WMG pipeline - `pytest -v tests/unit/wmg_processing`

### Run functional tests for WMG

Run functional tests for WMG api against the `dev` environment.
**NOTE**: `dev` environment is a remote environment. These functional tests run locally against a backend in a remote environment called `dev`.

1. `AWS_PROFILE=single-cell-dev DEPLOYMENT_STAGE=dev pytest -v tests/functional/backend/wmg/test_wmg_api.py`
