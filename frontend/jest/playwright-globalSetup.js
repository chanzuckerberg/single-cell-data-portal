const featureFlags = require("./featureFlags");
module.exports = async () => {
  process.env.STORAGE = JSON.stringify(featureFlags);
};
