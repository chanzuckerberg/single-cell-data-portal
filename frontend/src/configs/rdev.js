const {
  DEV_AUTH0_CLIENT_ID,
  PLAUSIBLE_DATA_DOMAIN_STAGING,
} = require("./common");

const configs = {
  API_URL: "$API_URL",
  AUTH0_CLIENT_ID: DEV_AUTH0_CLIENT_ID,
  AUTH0_DOMAIN: "login.cellxgene.dev.single-cell.czi.technology",
  PLAUSIBLE_DATA_DOMAIN: PLAUSIBLE_DATA_DOMAIN_STAGING,
};

if (typeof module !== "undefined") module.exports = configs;
