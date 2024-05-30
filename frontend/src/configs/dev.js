const { PLAUSIBLE_DATA_DOMAIN_STAGING } = require("./common");

const API_URL = "https://api.cellxgene.dev.single-cell.czi.technology";
const configs = {
  API_URL,
  DE_API_URL: API_URL,
  WMG_API_URL: API_URL,
  CELLGUIDE_DATA_URL:
    "https://cellguide.cellxgene.dev.single-cell.czi.technology",
  CENSUS_MODELS_DATA_URL: "https://contrib.cellxgene.cziscience.com",
  PLAUSIBLE_DATA_DOMAIN: PLAUSIBLE_DATA_DOMAIN_STAGING,
  SENTRY_DEPLOYMENT_ENVIRONMENT: "staging",
};

if (typeof module !== "undefined") module.exports = configs;
