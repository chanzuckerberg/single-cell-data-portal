# Chalice App Template

To deploy the app, type `make deploy` in this directory.

Filename                  | Purpose                           | Information links
--------------------------|-----------------------------------|------------------------------------------
`app.py`                  |The application entry point        | [Chalice Docs](https://chalice.readthedocs.io/en/latest/)
`requirements-dev.txt`    |Developer environment dependencies | [Pip requirements files](https://pip.readthedocs.io/en/1.1/requirements.html)
`requirements.txt`        |Application dependencies           | [Chalice App Packaging](https://chalice.readthedocs.io/en/latest/topics/packaging.html)
`Makefile`                |Tools for packaging and deploying  | [Automation and Make](https://swcarpentry.github.io/make-novice/)
`.chalice/config.in.json`	|A template for the Chalice config file for the app    | [Chalice Configuration File](https://chalice.readthedocs.io/en/latest/topics/configfile.html)
`iam/policy-template/browser-lamnbda.json`|IAM policy for the app's IAM role  | [Lambda Permissions](https://docs.aws.amazon.com/lambda/latest/dg/intro-permission-model.html)
`test/test.py`            |Test suite template                | [Python unittest](https://docs.python.org/3/library/unittest.html)

## How to create a new app from this template
1. Install the dependencies: `pip install -r requirements-dev.txt` and Terraform (`brew install terraform` or
   https://www.terraform.io/downloads.html)
1. Configure the AWS CLI (`pip install awscli`; `aws configure`).
1. `source ./environment` to set the correct environment variables.
1. Edit `.chalice/config.in.json` to set the name of your app and Lambda settings like memory, timeout, reserved
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
IAM role is set from the file `iam/policy-template/{app_name}-lambda.json`. Edit this policy and repeat the deployment
if your Lambda needs access to other AWS APIs. You can also edit the Makefile to parameterize this file or generate 
it from a template as needed. (The setting `autogen_policy` must be set to `false` in `.chalice/config.json` for 
Chalice to use this file.)
