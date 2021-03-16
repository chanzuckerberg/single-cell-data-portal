const {
  SecretsManagerClient,
  GetSecretValueCommand,
} = require("@aws-sdk/client-secrets-manager");
const featureFlags = require("./featureFlags");

const endpoint = new URL(process.env["BOTO_ENDPOINT_URL"]);

const client = new SecretsManagerClient({
  endpoint,
});

const secretValueRequest = {
  SecretId: "corpora/backend/dev/auth0-secret",
};

const command = new GetSecretValueCommand(secretValueRequest);

module.exports = async () => {
  try {
    console.log(await client.send(command));
    const secret = JSON.parse((await client.send(command)).SecretString);
    process.env.TEST_ACCOUNT_PASS = secret.test_account_password;
  } catch (error) {
    console.error(error);
  }
  process.env.STORAGE = JSON.stringify(featureFlags);
};
