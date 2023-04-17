const { PLAUSIBLE_DATA_DOMAIN_STAGING } = require("./common");

const configs = {
  API_URL: "https://api.cellxgene.dev.single-cell.czi.technology",
  PLAUSIBLE_DATA_DOMAIN: PLAUSIBLE_DATA_DOMAIN_STAGING,
  SENTRY_DEPLOYMENT_ENVIRONMENT: "staging",
};

if (typeof module !== "undefined") module.exports = configs;
