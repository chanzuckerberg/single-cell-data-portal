const { PLAUSIBLE_DATA_DOMAIN_STAGING } = require("./common");

const configs = {
  API_URL: "https://canonical-backend.rdev.single-cell.czi.technology",
  CELLGUIDE_DATA_URL:
    "https://cellguide.cellxgene.dev.single-cell.czi.technology",
  PLAUSIBLE_DATA_DOMAIN: PLAUSIBLE_DATA_DOMAIN_STAGING,
  SENTRY_DEPLOYMENT_ENVIRONMENT: "staging",
};

if (typeof module !== "undefined") module.exports = configs;
