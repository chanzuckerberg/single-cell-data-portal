
resource "aws_secretsmanager_secret" "auth0" {
  name = "dcp/backend/browser/test/auth0-secret"
  tags = {
    project   = var.project
    env       = var.env
    service   = var.service
    owner     = var.owner
    managedBy = "terraform"
  }
}


data "aws_iam_policy_document" "policy" {
  statement {
    sid    = "ReadAuth0TestSecret"
    effect = "Allow"
    actions = [
      "secretsmanager:GetSecretValue",
      "secretsmanager:DescribeSecret"
    ]

    resources = [
      aws_secretsmanager_secret.auth0.arn,
    ]
  }
}

resource "aws_iam_policy" "github_actions" {
  name        = "github_actions"
  path        = "/"
  policy      = data.aws_iam_policy_document.policy.json
  description = "Provides github actions access to aws for running tests."
}
