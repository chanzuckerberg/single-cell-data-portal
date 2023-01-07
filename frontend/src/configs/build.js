const { PLAUSIBLE_DATA_DOMAIN_STAGING } = require("./common");

const configs = {
  API_URL: "RUNTIME_VAR_API_URL",
  PLAUSIBLE_DATA_DOMAIN: "RUNTIME_VAR_PLAUSIBLE_DATA_DOMAIN",
  SENTRY_DEPLOYMENT_ENVIRONMENT: "RUNTIME_VAR_DEPLOYMENT_ENVIRONMENT",
};

if (typeof module !== "undefined") module.exports = configs;
