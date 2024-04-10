import os

from backend.common.utils.secret_config import SecretConfig


class CorporaConfig(SecretConfig):
    environ_source = "CORPORA_CONFIG"

    def __init__(self, *args, **kwargs):
        super().__init__("backend", secret_name="config", **kwargs)

    def get_defaults_template(self):

        # Set collections_base_url
        rdev_prefix = os.environ.get("REMOTE_DEV_PREFIX")
        dataset_assets_base_url = None
        if rdev_prefix:
            collections_base_url = f"https://{rdev_prefix.strip('/')}-frontend.rdev.single-cell.czi.technology"
            dataset_assets_base_url = f"https://env-rdev-datasets.s3.us-west-2.amazonaws.com/{rdev_prefix.strip('/')}"
        else:
            deployment_stage = os.environ.get("DEPLOYMENT_STAGE")
            if deployment_stage == "test":
                collections_base_url = "https://frontend.corporanet.local:3000"
            elif deployment_stage == "prod":
                collections_base_url = "https://cellxgene.cziscience.com"
            else:
                collections_base_url = f"https://cellxgene.{deployment_stage}.single-cell.czi.technology"

        template = {
            "upload_max_file_size_gb": 30,
            "submission_bucket": os.getenv("DATASET_SUBMISSIONS_BUCKET", "cellxgene-dataset-submissions-test"),
            "collections_base_url": collections_base_url,
            "dataset_assets_base_url": dataset_assets_base_url,
            "ingest_memory_modifier": 3,  # adds (x*100)% memory overhead
            "ingest_min_vcpu": 2,
            # The largest machine we are allocating is r5a.24xlarge. This machine has 768GB of memory and 96 vCPUs.
            "ingest_max_vcpu": 96,  # assume 8GB per vCPU
            "ingest_swap_modifier": 0,  # 0 b/c no swap machines are used.
            "ingest_max_swap_memory_mb": 300000,
        }
        upload_snf_arn = os.getenv("UPLOAD_SFN_ARN")
        if upload_snf_arn:
            template["upload_sfn_arn"] = upload_snf_arn
        return template


class CorporaDbConfig(SecretConfig):
    environ_source = "CORPORA_DB_CONFIG"

    def __init__(self, *args, **kwargs):
        super().__init__(
            component_name="backend",
            secret_name="database",
            **kwargs,
        )

    def get_defaults_template(self):
        # The db secret for remote dev envs is {"remote_dev_uri": "postgresql://blah"}
        # instead of {"database_uri": "postgresql://blah"} so we can add a suffix here
        # based on the remote dev env name.
        remote_dev_prefix = os.getenv("REMOTE_DEV_PREFIX", "")
        if not remote_dev_prefix:
            return {}
        return {
            "database_uri": "{remote_dev_uri}" + remote_dev_prefix,
        }


class CorporaAuthConfig(SecretConfig):
    """
    For a description of the secret key contents, see backend/config/auth0-secret-template.json.
    """

    environ_source = "CORPORA_AUTH_CONFIG"

    def __init__(self, *args, **kwargs):
        super().__init__(
            component_name="backend",
            secret_name="auth0-secret",
            **kwargs,
        )

    def get_defaults_template(self):
        template = {
            "api_authorize_url": "{api_base_url}/authorize",
            "api_token_url": "{api_base_url}/oauth/token",
            "api_userinfo_url": "{api_base_url}/userinfo",
            "internal_url": "{api_base_url}",
            "issuer": [],
            "retry_status_forcelist": [429, 500, 502, 503, 504],
        }
        template["issuer"].append(self.api_base_url + "/" if not self.api_base_url.endswith("/") else self.api_base_url)
        template["issuer"].append("https://" + self.auth0_domain + "/")
        if self.config.get("api_signin_url"):
            # Adding the API sign in URL to the list of allow token issues. This allow the API to accept Auth token
            # generated for testing. Used in dev and staging.
            template["issuer"].append(self.api_signin_url)

        return template


class CorporaCloudfrontConfig(SecretConfig):
    environ_source = "CORPORA_CLOUDFRONT_CONFIG"

    def __init__(self, *args, **kwargs):
        super().__init__("backend", secret_name="cloudfront", **kwargs)
