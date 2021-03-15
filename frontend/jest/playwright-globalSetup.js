const {
  SecretsManagerClient,
  GetSecretValueCommand,
} = require("@aws-sdk/client-secrets-manager");
const featureFlags = require("./featureFlags");

const client = new SecretsManagerClient({
  endpoint: process.env["BOTO_ENDPOINT_URL"],
});

const secretValueRequest = {
  SecretId: "corpora/backend/dev/auth0-secret['test_account_password']",
};

const command = new GetSecretValueCommand(secretValueRequest);

module.exports = async () => {
  try {
    const secret = await client.send(command);
    process.env.TEST_ACCOUNT_PASS = secret;
  } catch (error) {
    console.error(error);
  }
  process.env.STORAGE = JSON.stringify(featureFlags);
};
