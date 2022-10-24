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

/**
 * (thuang): From oauth/users.json
 * `frontend/tests/common/constants.ts` `API_URL` needs to be `"http://backend.corporanet.local:5000"`
 * If `API_URL` is dev, staging, prod, or rdev, use their corresponding OIDC credentials
 */
const LOCAL_CONTAINER_TEST_USERNAME = "User1";
const LOCAL_CONTAINER_TEST_PASSWORD = "pwd";

const DEPLOYED_TEST_USERNAME = "user@example.com";

/**
 * (thuang): `process.env.TEST_ACCOUNT_USER` and `process.env.TEST_ACCOUNT_PASS`
 * are exported from `frontend/playwright-globalSetup.ts`
 */
const deployedTestUsername =
  process.env.TEST_ACCOUNT_USER || DEPLOYED_TEST_USERNAME;

const USERNAMES = {
  dev: deployedTestUsername,
  happy: LOCAL_CONTAINER_TEST_USERNAME,
  local: LOCAL_CONTAINER_TEST_USERNAME,
  localProd: LOCAL_CONTAINER_TEST_USERNAME,
  prod: deployedTestUsername,
  rdev: deployedTestUsername,
  staging: deployedTestUsername,
};

export const TEST_USERNAME = USERNAMES[TEST_ENV] || deployedTestUsername;

export const TEST_PASSWORD =
  process.env.TEST_ACCOUNT_PASS || LOCAL_CONTAINER_TEST_PASSWORD || "";

export const TEST_URL = TEST_ENV_TO_TEST_URL[TEST_ENV];

export const BLUEPRINT_SAFE_TYPE_OPTIONS = { delay: 50 };

export const COOKIE_SESSION = process.env.COOKIE_SESSION || "";
export const COOKIE_CXG_USER = process.env.COOKIE_CXG_USER || "";
