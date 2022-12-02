const { PLAUSIBLE_DATA_DOMAIN_STAGING } = require("./common");

// (thuang): For local development, please copy the content of this file
// to a new file named `configs.js` in this directory.

const configs = {
  API_URL: "https://backend.corporanet.local:5000",
  AUTH0_URL: "https://oidc.corporanet.local:8443",
  PLAUSIBLE_DATA_DOMAIN: PLAUSIBLE_DATA_DOMAIN_STAGING,
};

if (typeof module !== "undefined") module.exports = configs;
