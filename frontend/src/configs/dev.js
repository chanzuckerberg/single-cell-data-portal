const { PLAUSIBLE_DATA_DOMAIN_STAGING } = require("./common");

const configs = {
  API_URL: "https://api.cellxgene.dev.single-cell.czi.technology",
  CELLGUIDE_DATA_URL:
    "https://cellguide.cellxgene.dev.single-cell.czi.technology",
  CENSUS_SPOTLIGHT_DATA_URL: "https://contrib.cellxgene.cziscience.com",
  PLAUSIBLE_DATA_DOMAIN: PLAUSIBLE_DATA_DOMAIN_STAGING,
  SENTRY_DEPLOYMENT_ENVIRONMENT: "staging",
};

if (typeof module !== "undefined") module.exports = configs;
