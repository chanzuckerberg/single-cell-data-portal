const TEST_URL = require("../tests/common/constants").TEST_URL;
const COOKIE_SESSION = require("../tests/common/constants").COOKIE_SESSION;
const COOKIE_CXG_USER = require("../tests/common/constants").COOKIE_CXG_USER;

const cookies = [];

if (COOKIE_SESSION) {
  cookies.push({
    domain: "localhost",
    expires: 1660161952,
    httpOnly: true,
    name: "session",
    path: "/",
    sameSite: "None",
    secure: true,
    session: true,
    size: 106,
    value: COOKIE_SESSION,
  });
}

if (COOKIE_CXG_USER) {
  cookies.push({
    domain: "api.cellxgene.staging.single-cell.czi.technology",
    expires: 1660161952,
    httpOnly: true,
    name: "cxguser",
    path: "/",
    sameSite: "None",
    secure: true,
    session: false,
    size: 1711,
    value: COOKIE_CXG_USER,
  });
}

module.exports = {
  cookies: [
    ...cookies,
    // (thuang) Uncomment and add your staging cookies here to skip login when
    // running E2E tests locally IF your local FE is connecting to staging API!
    // {
    //   domain: "api.cellxgene.staging.single-cell.czi.technology",
    //   expires: 1660161952,
    //   httpOnly: true,
    //   name: "cxguser",
    //   path: "/",
    //   sameSite: "None",
    //   secure: true,
    //   session: false,
    //   size: 1711,
    // // Add your staging cxguser cookie here
    //   value: "",
    // },
    // {
    //   domain: "localhost",
    //   expires: 1660161952,
    //   httpOnly: true,
    //   name: "session",
    //   path: "/",
    //   sameSite: "None",
    //   secure: true,
    //   session: true,
    //   size: 106,
    // // Add your staging session cookie here
    //   value: "",
    // },
  ],
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
