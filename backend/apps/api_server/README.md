# Portal Backend

This application serves as the API for the Data Portal website.

# App Overview

Filename                  | Purpose                           | Information links
--------------------------|-----------------------------------|------------------------------------------
`app.py`                  |The application entry point        |
`requirements.txt`        |Application dependencies           |
`requirements-dev.txt`    |Developer environment dependencies | [Pip requirements files](https://pip.readthedocs.io/en/1.1/requirements.html)
`Makefile`                |Tools for packaging and deploying  | [Automation and Make](https://swcarpentry.github.io/make-novice/)

## Development
See [DEV_ENV.md](../../../DEV_ENV.md) for a local dev quickstart.

## Deploy
1. Set `DEPLOYMENT_STAGE` to a new or existing stage: `test`, `staging`, `prod`
1. Install dependencies `pip install -r requirements-dev.txt`
1. Configure the AWS CLI
    1. `pip install awscli`
    1. `aws configure`
1. [Optional] Edit `app.py` and `requirements.txt` to modify the app.
1. Deploy your app:
   - For creating a new deployment run `make deploy`.
   - For updating an existing deployment run `make ci-deploy`.
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
it from a template as needed.
