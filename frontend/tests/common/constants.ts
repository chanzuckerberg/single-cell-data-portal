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
  happy: "http://frontend.corporanet.local:3000",
  local: "http://localhost:3000",
  localProd: "http://localhost:9000",
  prod: "https://cellxgene.cziscience.com",
  rdev: process.env.RDEV_LINK || "",
  dev: "https://cellxgene.dev.single-cell.czi.technology",
  staging: "https://cellxgene.staging.single-cell.czi.technology",
};

export const TEST_URL = TEST_ENV_TO_TEST_URL[TEST_ENV];
export const TEST_USERNAME =
  TEST_ENV === "happy"
    ? (process.env.TEST_ACCOUNT_USER as string)
    : "user@example.com";
export const TEST_PASSWORD = process.env.TEST_ACCOUNT_PASS || "";

export const BLUEPRINT_SAFE_TYPE_OPTIONS = { delay: 50 };
