// (thuang): For local development, please copy the content of this file
// to a new file named `configs.js` in this directory.

// Staging
const API_URL = "https://api.cellxgene.staging.single-cell.czi.technology";
// Prod
// const API_URL = "https://api.cellxgene.cziscience.com";
// Local container
// const API_URL = "https://backend.corporanet.local:5000";

const configs = {
  API_URL,
  DE_API_URL: API_URL,
  WMG_API_URL: API_URL,
  CELLGUIDE_DATA_URL:
    "https://cellguide.cellxgene.dev.single-cell.czi.technology",
  CENSUS_MODELS_DATA_URL: "https://contrib.cellxgene.cziscience.com",
  PLAUSIBLE_DATA_DOMAIN: "cellxgene.staging.single-cell.czi.technology",
};

if (typeof module !== "undefined") module.exports = configs;
