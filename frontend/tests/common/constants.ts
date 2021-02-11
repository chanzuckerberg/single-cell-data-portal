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
  switch (env) {
    case "dev":
    case "rdev":
      return "cellxgene-smoke-test+dev@chanzuckerberg.com";
    case "staging":
      return "cellxgene-smoke-test+staging@chanzuckerberg.com";
    case "prod":
      return "cellxgene-smoke-test+prod@chanzuckerberg.com";
    default:
      return process.env.AUTH_EMAIL || "dev";
  }
};

export const TEST_URL = TEST_ENV_TO_TEST_URL[TEST_ENV];
export const TEST_EMAIL = TEST_ENV_TO_TEST_EMAIL(TEST_ENV);
