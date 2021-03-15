const {
  SecretsManagerClient,
  GetSecretValueCommand,
} = require("@aws-sdk/client-secrets-manager");
const featureFlags = require("./featureFlags");

const client = new SecretsManagerClient({
  endpoint: process.env["BOTO_ENDPOINT_URL"],
  region: "us-west-2",
});

const secretValueRequest = {
  SecretId: "corpora/backend/dev/auth0-secret",
};

const command = new GetSecretValueCommand(secretValueRequest);

module.exports = async () => {
  try {
    const secret = JSON.parse((await client.send(command)).SecretString);
    process.env.TEST_ACCOUNT_PASS = secret.test_account_password;
  } catch (error) {
    console.error(error);
  }
  process.env.STORAGE = JSON.stringify(featureFlags);
};
