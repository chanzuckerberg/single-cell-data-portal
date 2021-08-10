const TEST_URL = require("../tests/common/constants").TEST_URL;

module.exports = {
  cookies: [],
  origins: [
    {
      localStorage: [
        { name: "cxg-ff-gs", value: "yes" },
        { name: "cxg-ff-curator", value: "yes" },
        { name: "cxg-ff-rc", value: "yes" },
      ],
      origin: TEST_URL,
    },
  ],
};
