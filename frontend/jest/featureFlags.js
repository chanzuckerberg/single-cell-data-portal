const TEST_URL = require("../tests/common/constants").TEST_URL;

module.exports = {
  cookies: [],
  origins: [
    {
      localStorage: [
        { name: "cxg-ff-auth", value: "yes" },
        { name: "cxg-ff-cc", value: "yes" },
      ],
      origin: TEST_URL,
    },
  ],
};
