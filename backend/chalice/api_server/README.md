# Portal Backend

This application serves as the API for the Data Portal website.

# Chalice App Overview

Filename                  | Purpose                           | Information links
--------------------------|-----------------------------------|------------------------------------------
`./chalice/app.py`                  |The application entry point        | [Chalice Docs](https://chalice.readthedocs.io/en/latest/)
`./chalice/requirements.txt`        |Application dependencies           | [Chalice App Packaging](https://chalice.readthedocs.io/en/latest/topics/packaging.html)
`./chalice/.chalice/config.in.json`	|A template for the Chalice config file for the app    | [Chalice Configuration File](https://chalice.readthedocs.io/en/latest/topics/configfile.html)
`requirements-dev.txt`    |Developer environment dependencies | [Pip requirements files](https://pip.readthedocs.io/en/1.1/requirements.html)
`Makefile`                |Tools for packaging and deploying  | [Automation and Make](https://swcarpentry.github.io/make-novice/)
`../../config/iam-policy-template/corpora-api-lambda.json`|IAM policy for the app's IAM role  | [Lambda Permissions](https://docs.aws.amazon.com/lambda/latest/dg/intro-permission-model.html)

## Development
1.  Ensure project [pre-requisites](../../../README.md#Pre-requisites) are met.

1.  Set the following environment variables:

    ```shell
    export DEPLOYMENT_STAGE=dev
    export AWS_PROFILE=single-cell-dev
    ```

1.  Install dependencies

     ```shell
     pip install -r requirements.txt
     ```

1.  Enable local connection to the private `dev` RDS instance:

     ```shell
     ssh -L 5432:corpora-dev-corpora-api.cluster-c81u9imopfwl.us-west-2.rds.amazonaws.com:5432 bastion.dev.single-cell.czi.technology
     ```

    This command opens an ssh tunnel from `localhost:5432` to the RDS connection endpoint via the `bastion` server.
    The local port `5432` is fixed and encoded in the DB connection string stored in the AWS Secret at
    `corpora/backend/<DEPLOYMENT_STAGE>/database_local`.

    Note: Since the default PostgreSQL port is `5432`, the above command will conflict with a local PostgreSQL instance
    if one is runinng. The currently running instance will need to be stopped before successfully opening an ssh tunnel.

1.  Deploy the Chalice app to http://localhost:5000

     ```shell
     make local-server
     ```

## Deploy
1. Set `DEPLOYMENT_STAGE` to a new or existing stage: `dev`, `staging`, `prod`
1. Install dependencies `pip install -r requirements-dev.txt`
1. Configure the AWS CLI
    1. `pip install awscli`
    1. `aws configure`
1. [Optional] Edit `./chalice/.chalice/config.in.json` to set the name of your app and Lambda settings like memory, timeout, reserved
   concurrency, tags, and environment variables.
1. [Optional] Edit `./chalice/app.py` and `./chalice/requirements.txt` to modify the app.
1. Deploy your app:
   - For creating a new deployment run `make new-deploy`.
   - For updating an existing deployment run `make deploy`.
1. The deployment results, including your API's Endpoint URL, will be printed to the terminal.
1. If needed, assign
   a [Custom Domain Name](https://docs.aws.amazon.com/apigateway/latest/developerguide/how-to-custom-domains.html),
   [ACM certificate](https://aws.amazon.com/certificate-manager/) and [Route 53 CNAME](https://aws.amazon.com/route53/)
   to your app in the [API Gateway AWS Console](https://console.aws.amazon.com/apigateway/home#/custom-domain-names),
   so users can reach your app on a friendly domain like `https://app.czi.technology`.

To redeploy your app after updating, run `make deploy` again. To undeploy the app and delete all associated resources,
run `make destroy`.

## Testing
Tests are run in the top level directory `corpora-data-portal`.

See the [top level README](https://github.com/chanzuckerberg/corpora-data-portal/blob/main/README.md#commands)
for how to run tests.

## Managing the Lambda IAM role and assume role policy
Your Lambda function is assigned an IAM role that controls the permissions given to the Lambda's AWS credentials. This
IAM role is set from the file [corpora-api-lambda.json](../../config/iam-policy-templates/corpora-api-lambda.json). Edit this policy and repeat the deployment
if your Lambda needs access to other AWS APIs. You can also edit the Makefile to parameterize this file or generate 
it from a template as needed. (The setting `autogen_policy` must be set to `false` in `.chalice/config.json` for 
Chalice to use this file.)
