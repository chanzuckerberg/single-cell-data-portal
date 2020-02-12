#!/usr/bin/env -S jq -f

# This jq filter pipeline takes the Terraform script output of `chalice package --pkg-format terraform` and injects a
# few workarounds for issues that we found. It is invoked by the `make deploy` target.


# Instead of referencing a local Lambda zipfile and uploading it in the CreateFunction API call, upload it to S3 first.
# This helps deploy multiple Lambdas from the same distribution via slow network connections (the zipfile only gets
# uploaded once). Also, we found that uploading the zipfile in CreateFunction has a very non-zero API error rate.

del(.resource.aws_lambda_function[].filename)
| .resource.aws_lambda_function[].s3_bucket=env.TF_S3_BUCKET
| .resource.aws_lambda_function[].s3_key=env.LAMBDA_SHA+".zip"
| .resource.aws_lambda_function[].source_code_hash=env.LAMBDA_SHA


# Workaround for https://github.com/terraform-providers/terraform-provider-aws/issues/420.
# When updating an existing deployment while using an API Gateway custom domain name with a base path mapping, Terraform
# is unable to order the operations correctly and stops with the following error:
# "BadRequestException: Active stages pointing to this deployment must be moved or deleted".

| .resource.aws_api_gateway_deployment.rest_api.lifecycle.create_before_destroy=true


# Workaround for a problem described in https://forums.aws.amazon.com/thread.jspa?threadID=256140.

# AWS API Gateway is overly opinionated on how to handle content type negotiation. It also has a number of bugs in code
# related to this content type negotiation. However, API Gateway is the only suitable front-end for Lambda over HTTP
# (ALB has even worse limits on request/response sizes, and is not supported by Chalice), so we have to use it.

# By default, API Gateway will refuse to handle any binary responses to requests that don't list the specific response
# content type in the Accept header. This behavior can be fixed by using the "proxy integration" and adding "*/*" to
# the list of binary content types that API Gateway should leave alone. We implement this workaround in
# dcpquery.api.ChaliceWithGzipBinaryResponses
# (https://github.com/HumanCellAtlas/query-service/blob/fb0a533/dcpquery/api/__init__.py#L111).

# However, this workaround has an unexpected side effect of breaking the "mock integration" for OPTIONS requests, which
# are needed to support CORS preflight by browsers using POST/PUT/PATCH to talk to the application. "Mock integration"
# here means API Gateway never calls the Lambda, instead returning a static response that we specify. This is faster
# than calling the Lambda, and is supported by Chalice
# (https://chalice.readthedocs.io/en/latest/quickstart.html#tutorial-cors-support).
# But when we force a "binary" content type, API Gateway becomes unable to construct this static response and crashes.
# The workaround is to set the "mock integration" contentHandling property to "CONVERT_TO_TEXT" as described in the
# forum post linked above.

| .data.template_file.chalice_api_swagger.template|=(fromjson|(.paths[].options["x-amazon-apigateway-integration"]|select(.type=="mock")|.contentHandling)="CONVERT_TO_TEXT"|tojson)
