type TEST_ENV = "local" | "localProd" | "dev" | "staging" | "prod" | "rdev";

export const TEST_ENV: TEST_ENV = (process.env.TEST_ENV as TEST_ENV) || "dev";

const TEST_ENV_TO_TEST_URL = {
  dev: "https://cellxgene.dev.single-cell.czi.technology",
  local: "http://localhost:8000",
  localProd: "http://localhost:9000",
  prod: "https://cellxgene.cziscience.com",
  rdev: process.env.RDEV_LINK || "",
  staging: "https://cellxgene.staging.single-cell.czi.technology",
};

const TEST_ENV_TO_TEST_EMAIL = (env: TEST_ENV) => {
  return "";
};

export const TEST_URL = TEST_ENV_TO_TEST_URL[TEST_ENV];
export const TEST_EMAIL = TEST_ENV_TO_TEST_EMAIL(TEST_ENV);

export const BLUEPRINT_SAFE_TYPE_OPTIONS = { delay: 50 };
