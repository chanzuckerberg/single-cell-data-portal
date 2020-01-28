# Consumer Prototype

Serve artifacts to users.


### Install Dependencies

The consumer requires Python 3.7 to run.

Clone the repo and install dependencies:
```
git clona git@github.com:chanzuckerberg/dcp-prototype.git
cd dcp-prototype/
pip install -r requirements-dev.txt
cd dcp-prototype/dcp-prototype/consumer
pip install -r requirements.txt
```

Also install [terraform from Hashicorp](https://www.terraform.io/) from your favourite package manager.

# deploy infra

```bash
export AWS_PROFILE=czi-hca-dev
export AWS_DEFAULT_REGION=us-east-1
export AWS_DEFAULT_OUTPUT=json
cd dcp-prototype/
source ./dcp-prototype/environment
make -C infra plan COMPONENT=consumer
```