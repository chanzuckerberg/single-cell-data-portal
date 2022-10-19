type TEST_ENV =
  | "local"
  | "localProd"
  | "staging"
  | "prod"
  | "rdev"
  | "happy"
  | "dev";

export const TEST_ENV: TEST_ENV = (process.env.TEST_ENV as TEST_ENV) || "local";

const TEST_ENV_TO_TEST_URL = {
  dev: "https://cellxgene.dev.single-cell.czi.technology",
  happy: "http://frontend.corporanet.local:3000",
  local: "http://localhost:3000",
  localProd: "http://localhost:9000",
  prod: "https://cellxgene.cziscience.com",
  rdev: process.env.RDEV_LINK || "",
  staging: "https://cellxgene.staging.single-cell.czi.technology",
};

export const TEST_URL = TEST_ENV_TO_TEST_URL[TEST_ENV];

// (thuang): From oauth/users.json
const LOCAL_TEST_USERNAME = "User1";
const LOCAL_TEST_PASSWORD = "pwd";

const DEPLOYED_TEST_USERNAME = "user@example.com";

const USERNAMES = {
  dev: DEPLOYED_TEST_USERNAME,
  happy: LOCAL_TEST_USERNAME,
  local: LOCAL_TEST_USERNAME,
  localProd: LOCAL_TEST_USERNAME,
  prod: DEPLOYED_TEST_USERNAME,
  rdev: DEPLOYED_TEST_USERNAME,
  staging: DEPLOYED_TEST_USERNAME,
};

export const TEST_USERNAME = USERNAMES[TEST_ENV] || DEPLOYED_TEST_USERNAME;

export const TEST_PASSWORD =
  process.env.TEST_ACCOUNT_PASS || LOCAL_TEST_PASSWORD;

export const BLUEPRINT_SAFE_TYPE_OPTIONS = { delay: 50 };

export const COOKIE_SESSION = process.env.COOKIE_SESSION || "";
export const COOKIE_CXG_USER = process.env.COOKIE_CXG_USER || "";
