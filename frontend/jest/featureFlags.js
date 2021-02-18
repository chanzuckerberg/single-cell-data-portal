const TEST_URL = require("../tests/common/constants").TEST_URL;

module.exports = {
  cookies: [],
  origins: [
    {
      localStorage: [
        { name: "cxg-ff-auth", value: "true" },
        { name: "cxg-ff-cc", value: "true" },
      ],
      origin: TEST_URL,
    },
  ],
};
