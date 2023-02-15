// (thuang): For local development, please copy the content of this file
// to a new file named `configs.js` in this directory.

const configs = {
  API_URL: "https://backend.corporanet.local:5000",
  PLAUSIBLE_DATA_DOMAIN: "cellxgene.staging.single-cell.czi.technology"
};

if (typeof module !== "undefined") module.exports = configs;
