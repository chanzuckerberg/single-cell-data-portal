import {
  GetSecretValueCommand,
  SecretsManagerClient,
} from "@aws-sdk/client-secrets-manager";

let endpoint;

try {
  endpoint = new URL(process.env["BOTO_ENDPOINT_URL"] as string);
} catch (e) {
  console.log("BOTO_ENDPOINT_URL not assigned, assuming running on deployment");
  endpoint = "";
}

const client = new SecretsManagerClient({
  endpoint: String(endpoint),
});

const deployment_stage = process.env.DEPLOYMENT_STAGE || "test";

const secretValueRequest = {
  SecretId: `corpora/backend/${deployment_stage}/auth0-secret`,
};

const command = new GetSecretValueCommand(secretValueRequest);

module.exports = async () => {
  try {
    const secret = JSON.parse(
      (await client.send(command)).SecretString || "null"
    );
    process.env.TEST_ACCOUNT_USER = secret.test_account_username;
    process.env.TEST_ACCOUNT_PASS = secret.test_account_password;
  } catch (error) {
    console.error(error);
  }
};
