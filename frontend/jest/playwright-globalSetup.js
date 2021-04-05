const {
  SecretsManagerClient,
  GetSecretValueCommand,
} = require("@aws-sdk/client-secrets-manager");
const featureFlags = require("./featureFlags");
let endpoint;
try {
  endpoint = new URL(process.env["BOTO_ENDPOINT_URL"]);
} catch (e) {
  console.log("BOTO_ENDPOINT_URL not assigned, assuming running on deployment");
  endpoint = "";
}

const client = new SecretsManagerClient({
  endpoint,
});

const secretValueRequest = {
  SecretId: "corpora/backend/dev/auth0-secret",
};

const command = new GetSecretValueCommand(secretValueRequest);

module.exports = async () => {
  try {
    const secret = JSON.parse((await client.send(command)).SecretString);
    process.env.TEST_ACCOUNT_USER = secret.test_account_username;
    process.env.TEST_ACCOUNT_PASS = secret.test_account_password;
  } catch (error) {
    console.error(error);
  }
  process.env.STORAGE = JSON.stringify(featureFlags);
};
