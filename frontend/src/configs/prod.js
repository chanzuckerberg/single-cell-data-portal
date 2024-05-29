const configs = {
  API_URL: "https://api.cellxgene.cziscience.com",
  DE_API_URL: "https://de-api.cellxgene.cziscience.com",
  WMG_API_URL: "https://wmg-api.cellxgene.cziscience.com",
  CELLGUIDE_DATA_URL: "https://cellguide.cellxgene.cziscience.com",
  CENSUS_MODELS_DATA_URL: "https://contrib.cellxgene.cziscience.com",
  PLAUSIBLE_DATA_DOMAIN: "cellxgene.cziscience.com",
  SENTRY_DEPLOYMENT_ENVIRONMENT: "production",
};

if (typeof module !== "undefined") module.exports = configs;
