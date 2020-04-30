# Cloudfront-Invalidator

A chalice application used to invalidate the cache of a cloudfront distribution for
an s3 bucket hosting a website. When a change is detected in the s3 bucket, the
AWS Cloudfront distribution for that bucket is invalidated. This allowschanges
made to a website to become immediately visible.

# Chalice App Overview

Filename                  | Purpose                           | Information links
--------------------------|-----------------------------------|------------------------------------------
[./app.py](app.py)                  |The application entry point        | [Chalice Docs](https://chalice.readthedocs.io/en/latest/)
[./requirements.txt](requirements.txt)        |Application dependencies           | [Chalice App Packaging](https://chalice.readthedocs.io/en/latest/topics/packaging.html)
[./.chalice/config.in.json](.chalice/config.in.json)	|A template for the Chalice config file for the app    | [Chalice Configuration File](https://chalice.readthedocs.io/en/latest/topics/configfile.html)
[Makefile](Makefile)                |Tools for packaging and deploying  | [Automation and Make](https://swcarpentry.github.io/make-novice/)
[iam/policy-templates/cloudfront-invalidator-lambda.json](../../iam/policy-templates/cloudfront-invalidator-lambda.json)|IAM policy for the app's IAM role  | [Lambda Permissions](https://docs.aws.amazon.com/lambda/latest/dg/intro-permission-model.html)

## Development
1.  Ensure your `awscli` is configured with the
    [required credentials and profiles](https://github.com/chanzuckerberg/dcp-prototype#configuration).

1.  Set the following environment variables:

    ```shell
    export DEPLOYMENT_STAGE=dev
    export AWS_PROFILE=single-cell-dev
    ```

1.  Install dependencies

     ```shell
     pip install -r ../../../requirements-dev.txt
     pip install -r requirement.txt
     ```

## How to deploy and manage the Chalice app
1. Set environment variables:
	```shell
	export DEPLOYMENT_STAGE=dev
	export AWS_PROFILE=single-cell-dev
	```
	
1. Install dependencies:
     ```shell
     pip install -r ../../../requirements-dev.txt
     pip install -r requirement.txt
     pip install awscli
     ```
1. Configure the AWS CLI
    1. `aws configure`
1. [Optional] Edit [./.chalice/config.in.json](.chalice/config.in.json) to set the name of your app and Lambda settings like memory, timeout, reserved
   concurrency, tags, and environment variables.
1. [Optional] Edit [./app.py](app.py)   and [./requirements.txt](requirements.txt) to modify the app.
1. Deploy your app by running `make deploy`.

To redeploy your app after updating, run `make deploy` again. To undeploy the app and delete all associated resources,
run `make destroy`.

## Testing
Tests are run in the top level directory [dcp_prototpye](../../../README.md).

See the [top level README](../../../README.md#testing)
for how to run tests.

## Managing the Lambda IAM role and assume role policy
Your Lambda function is assigned an IAM role that controls the permissions given to the Lambda's AWS credentials. This
IAM role is set from the file iam/policy-templates/cloudfront-invalidator-lambda.json](../../iam/policy-templates/cloudfront-invalidator-lambda.json). Edit this policy and repeat the deployment
if your Lambda needs access to other AWS APIs. You can also edit the Makefile to parameterize this file or generate 
it from a template as needed. (The setting `autogen_policy` must be set to `false` in `.chalice/config.json` for 
Chalice to use this file.)
