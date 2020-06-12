# Corpora: AWS CLI Configuration

Follow these steps to install and configure the `awscli`:
1. `pip install awscli`
1. Run `aws configure` or manually configure `czi-id` credentials:

    ```shell
    # ~/.aws/credentials example

    [czi-id]
    aws_access_key_id = ACCESS_KEY_ID
    aws_secret_access_key = SECRET_ACCESS_KEY
    ```

1.  [Configure](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-files.html)
    the `single-cell-dev` and/or `single-cell-prod` AWS profiles.

    ```shell
    # ~/.aws/config example

    [profile single-cell-dev]
    role_arn = arn:aws:iam::ACCOUNT_ID:role/poweruser
    source_profile = czi-id
    region = us-east-1

    [profile single-cell-prod]
    role_arn = arn:aws:iam::ACCOUNT_ID:role/poweruser
    source_profile = czi-id
    region = us-east-1
    ```

    Please contact #help-infra on Slack if you require access to `single-cell-dev` or `single-cell-prod`.
