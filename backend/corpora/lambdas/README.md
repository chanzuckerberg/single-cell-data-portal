# Lambdas

This directory contain code for AWS Lambdas deployed using [chalice](https://chalice.readthedocs.io/en/latest/).

# Deploy
1. Setup environment variables:
   * `DEPLOYMENT_STAGE`= the deployment to target
   * `AWS_PROFILE`= the AWS account to target
1. Run `make deploy APP_NAME=${app_name}` where `app_name` is the name of an application the [./chalice](chalice)

# Adding an Application
1. Run `chalice create-project` in [./chalice](chalice)
1. Create your app.
1. Add a makefile with the targets `deploy`, `destroy`, and `clean`

Notes:
* As much of the code should be imported from [./code](code) to allow to simplify unittesting.
* Store custom IAM policies in [./iam-policy-templates](iam-policy-templates/policy-templated]) with the name 
  `${APP_NAME}-lambda.json`.
