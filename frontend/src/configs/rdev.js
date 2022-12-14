const { PLAUSIBLE_DATA_DOMAIN_STAGING } = require("./common");

const configs = {
  API_URL: "$API_URL",
  AUTH0_CLIENT_ID: "R6GwxkCT5JwI8udLdlOtL3reVvy41QW5",
  AUTH0_DOMAIN: "login.cellxgene.dev.single-cell.czi.technology",
  PLAUSIBLE_DATA_DOMAIN: PLAUSIBLE_DATA_DOMAIN_STAGING,
};

if (typeof module !== "undefined") module.exports = configs;
