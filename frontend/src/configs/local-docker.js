const { PLAUSIBLE_DATA_DOMAIN_STAGING } = require("./common");

// Temporarily route API calls to the backend app to http://backend
// instead of http://localhost
const configs = {
  API_URL: "http://backend:5005",
  PLAUSIBLE_DATA_DOMAIN: PLAUSIBLE_DATA_DOMAIN_STAGING,
};

if (typeof module !== "undefined") module.exports = configs;
