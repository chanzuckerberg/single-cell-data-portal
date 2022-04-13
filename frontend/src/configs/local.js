const { PLAUSIBLE_DATA_DOMAIN_STAGING } = require("./common");

// (thuang): For local development, please copy the content of this file
// to a new file named `configs.js` in this directory.

const configs = {
  API_URL: "http://backend.corporanet.local:5000",
  PLAUSIBLE_DATA_DOMAIN: PLAUSIBLE_DATA_DOMAIN_STAGING,
};

if (typeof module !== "undefined") module.exports = configs;
