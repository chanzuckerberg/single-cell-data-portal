// (thuang): For local development, please copy the content of this file
// to a new file named `configs.js` in this directory.

const configs = {
  // Dev
  API_URL: "https://api.cellxgene.dev.single-cell.czi.technology",
  CELLGUIDE_DATA_URL:
    "https://cellguide.cellxgene.dev.single-cell.czi.technology",
  // Staging
  // API_URL: "https://api.cellxgene.staging.single-cell.czi.technology",
  // Prod
  // API_URL: "https://api.cellxgene.cziscience.com",
  // Local container
  // API_URL: "https://backend.corporanet.local:5000",
  PLAUSIBLE_DATA_DOMAIN: "cellxgene.staging.single-cell.czi.technology",
};

if (typeof module !== "undefined") module.exports = configs;
