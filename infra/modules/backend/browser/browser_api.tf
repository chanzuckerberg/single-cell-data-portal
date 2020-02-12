resource "aws_api_gateway_rest_api" "browser_backend" {
  name        = "${var.deployment_stage}-dcp-browser-backend"
  description = "This API is for serving the DCP 2.0 Browser Frontend"
}

resource "aws_api_gateway_resource" "browser_backend_" {
  rest_api_id = "${aws_api_gateway_rest_api.MyDemoAPI.id}"
  parent_id   = "${aws_api_gateway_rest_api.MyDemoAPI.root_resource_id}"
  path_part   = ""
}
