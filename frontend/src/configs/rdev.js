const { PLAUSIBLE_DATA_DOMAIN_STAGING } = require("./common");

const configs = {
  API_URL: "$API_URL",
  PLAUSIBLE_DATA_DOMAIN: PLAUSIBLE_DATA_DOMAIN_STAGING,
};

if (typeof module !== "undefined") module.exports = configs;
