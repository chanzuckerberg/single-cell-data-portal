const configs = {
  API_URL: "$API_URL",
  // plausible doesn't expect  a protocol in the domain
  PLAUSIBLE_DATA_DOMAIN: "$FRONTEND_URL".replace(/^https?:\/\//, ""),
};

if (typeof module !== "undefined") module.exports = configs;
