import { COOKIE_CXG_USER, COOKIE_SESSION, TEST_URL } from "./constants";

const cookies = [];

type SameSite = "Strict" | "Lax" | "None";

if (COOKIE_SESSION) {
  cookies.push({
    domain: "localhost",
    expires: 1660161952,
    httpOnly: true,
    name: "session",
    path: "/",
    sameSite: "None" as SameSite,
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
    sameSite: "None" as SameSite,
    secure: true,
    session: false,
    size: 1711,
    value: COOKIE_CXG_USER,
  });
}

const featureFlags = {
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
    //   sameSite: "None" as SameSite,
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
    //   sameSite: "None" as SameSite,
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

export default featureFlags;
