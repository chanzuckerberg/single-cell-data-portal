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
  happy: "https://frontend.corporanet.local:3000",
  local: "https://localhost:3000",
  localProd: "http://localhost:9000",
  prod: "https://cellxgene.cziscience.com",
  rdev: process.env.RDEV_LINK || "",
  staging: "https://cellxgene.staging.single-cell.czi.technology",
};

export const TEST_URL = TEST_ENV_TO_TEST_URL[TEST_ENV];

export const BLUEPRINT_SAFE_TYPE_OPTIONS = { delay: 50 };

export const SKIP_LOGIN = process.env.SKIP_LOGIN === "true" || false;

export const LOGIN_STATE_FILENAME = "loginState.json";
export const SHARED_LINK_FILTER = `${TEST_URL}/gene-expression?compare=disease&sexes=PATO%3A0000383&tissues=blood%2Clung&genes=DPM1%2CTNMD%2CTSPAN6&ver=2`;