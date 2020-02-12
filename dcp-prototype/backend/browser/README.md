# Chalice App Template

This repo is an AWS Lambda serverless app template with Terraform deployment management. It is intended to guide you to
good choices for ease of development and operations, without being overly opinionated about your workflow.

[AWS Lambda](https://aws.amazon.com/lambda/) is a serverless application platform that we use to build web applications,
services/APIs, and event handlers on AWS. Lambda scales quickly to accommodate burst loads, and integrates tightly
with [many AWS event sources](http://docs.aws.amazon.com/lambda/latest/dg/invoking-lambda-function.html), allowing us to
build simpler, more efficient and maintainable applications on AWS.

[Terraform](https://www.terraform.io/) is
an [Infrastructure-as-Code](https://en.wikipedia.org/wiki/Infrastructure_as_Code) framework that we use to manage our
cloud resources. Terraform allows us to automatically, reproducibly, and securely manage our AWS infrastructure.

Out of the box, Terraform and Lambda require a lot of configuration and domain knowledge to work together. You need to
know how to package your function and its dependencies, configure its IAM role, set up event sources, configure API
Gateway, set up a Terraform credentials source and state backend, and more, depending on your application.

In this repository, we provide the tools to take care of all this in a minimalistic, well-documented, intuitive way, so
you can focus on building your app.

To manage the Lambda packaging, we use [Chalice](https://github.com/aws/chalice), a Python microframework developed by
AWS for Lambda. The template in this repo builds upon the output of `chalice new-project`.

Instead of letting Chalice directly deploy our app, we run `chalice package`, which builds the app and produces
a Terraform module, `chalice.tf.json`. We then use Terraform to deploy the app from this module.

We also provide a recipe to easily manage Terraform credentials and state on AWS. The AWS API credentials are imported
on demand from your [AWS CLI](https://aws.amazon.com/cli/) config, which allows you to centrally manage your
credentials, [Assume Role](https://docs.aws.amazon.com/cli/latest/userguide/cli-roles.html) configurations, and regions
using standard AWS conventions. Terraform state files are saved using
the [Terraform S3 backend](https://www.terraform.io/docs/backends/types/s3.html), keyed by your app name.

To deploy the app, type `make deploy` in this directory.

Filename                  | Purpose                           | Information links
--------------------------|-----------------------------------|------------------------------------------
`app.py`                  |The application entry point        | [Chalice Docs](https://chalice.readthedocs.io/en/latest/)
`requirements-dev.txt`    |Developer environment dependencies | [Pip requirements files](https://pip.readthedocs.io/en/1.1/requirements.html)
`requirements.txt`        |Application dependencies           | [Chalice App Packaging](https://chalice.readthedocs.io/en/latest/topics/packaging.html)
`Makefile`                |Tools for packaging and deploying  | [Automation and Make](https://swcarpentry.github.io/make-novice/)
`.chalice/config.json`    |Chalice config file for the app    | [Chalice Configuration File](https://chalice.readthedocs.io/en/latest/topics/configfile.html)
`.chalice/policy-dev.json`|IAM policy for the app's IAM role  | [Lambda Permissions](https://docs.aws.amazon.com/lambda/latest/dg/intro-permission-model.html)
`test/test.py`            |Test suite template                | [Python unittest](https://docs.python.org/3/library/unittest.html)
`.travis.yml`             |Travis CI (CI/CD) configuration    | [Travis CI](https://docs.travis-ci.com/user/customizing-the-build/)

## How to create a new app from this template
1. Install the dependencies: `pip install -r requirements-dev.txt` and Terraform (`brew install terraform` or
   https://www.terraform.io/downloads.html)
1. Configure the AWS CLI (`pip install awscli`; `aws configure`).
1. Ensure the S3 bucket `tfstate-<YOUR_AWS_ACCOUNT_ID>` exists, or modify the Makefile to reference a different bucket.
1. Fork or copy the contents of this repo to a new directory.
1. Edit `.chalice/config.json` to set the name of your app and Lambda settings like memory, timeout, reserved
   concurrency, tags, and environment variables.
1. Edit `app.py` and `requirements.txt` to create your app.
1. Deploy your app by running `make deploy`. The deployment results, including your Lambda's EndpointURL, will be
   printed to the terminal. You can immediately test your app by running (for example)
   `http https://your-api-id.execute-api.us-east-1.amazonaws.com/api/` or opening the EndpointURL in a browser.
1. If needed, assign
   a [Custom Domain Name](https://docs.aws.amazon.com/apigateway/latest/developerguide/how-to-custom-domains.html),
   [ACM certificate](https://aws.amazon.com/certificate-manager/) and [Route 53 CNAME](https://aws.amazon.com/route53/)
   to your app in the [API Gateway AWS Console](https://console.aws.amazon.com/apigateway/home#/custom-domain-names),
   so users can reach your app on a friendly domain like `https://app.czi.technology`.

To redeploy your app after updating, run `make deploy` again. To undeploy the app and delete all associated resources,
run `make destroy`.

## Testing
The test suite in `test/test.py` runs Chalice in local mode for unit testing. You can invoke it using `make test`. This
test is also configured to run on [Travis CI](https://travis-ci.com).

## Managing the Lambda IAM role and assume role policy
Your Lambda function is assigned an IAM role that controls the permissions given to the Lambda's AWS credentials. This
IAM role is set from the file `.chalice/policy-dev.json`. Edit this policy and repeat the deployment if your Lambda
needs access to other AWS APIs. You can also edit the Makefile to parameterize this file or generate it from a template
as needed. (The setting `autogen_policy` must be set to `false` in `.chalice/config.json` for Chalice to use this file.)

## Handling other AWS event sources
See [Lambda Event Sources](https://chalice.readthedocs.io/en/latest/topics/events.html) in the Chalice docs for details
on how to connect your app to
[other AWS event sources](http://docs.aws.amazon.com/lambda/latest/dg/invoking-lambda-function.html). If the built-in
event handlers are not sufficient, you could also consider using the
[Domovoi project](https://github.com/kislyuk/domovoi), which builds upon Chalice to handle many other event source
configurations (but is not currently compatible with the Terraform workflow shown here).

## Monitoring your app
Lambda is automatically set up to emit logs to
[CloudWatch Logs](https://docs.aws.amazon.com/AmazonCloudWatch/latest/logs/WhatIsCloudWatchLogs.html), which you can
browse in the [AWS console](https://console.aws.amazon.com/cloudwatch/home#logs:) by selecting the log group for your
Lambda. You can see built-in CloudWatch metrics (invocations, errors, etc.) by selecting your Lambda in the
[Lambda AWS console](https://console.aws.amazon.com/lambda/home#/functions) and going to the **Monitoring** tab.
To tail and filter logs on the command line, you can use the `logs` command in the
[Aegea](https://github.com/kislyuk/aegea) package, for example:

    aegea logs /aws/lambda/my-lambda-APIHandler-ABCDEXAMPLE --start-time=-15m

[![Build Status](https://travis-ci.com/chanzuckerberg/chalice-app-template.svg?token=iPJHxi7MxMYqJkBxfGCC&branch=master)](https://travis-ci.com/chanzuckerberg/chalice-app-template)


##
- brew install jq moreutils gettext
