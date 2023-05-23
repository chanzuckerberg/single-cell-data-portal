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

export const ADD_TISSUE_BTN = "add-tissue-btn";
export const ADD_TISSUE_LBL = "get-started-step-1";
export const ADD_GENE_BTN = "add-gene-btn";
export const ADD_GENE_LBL = "get-started-step-2";
export const HOMO_SAPIENS_TERM_ID = "NCBITaxon:9606";
export const GENE_LABELS_ID = "[data-testid^=gene-label-]";
export const CELL_TYPE_LABELS_ID = "cell-type-name";
export const ADD_TISSUE_ID = "add-tissue-btn";
export const ADD_GENE_ID = "add-gene-btn";
export const GENE_DELETE_BUTTON = "gene-delete-button";
export const SOURCE_DATA_BUTTON_ID = "source-data-button";
export const SOURCE_DATA_LIST_SELECTOR = `[data-testid="source-data-list"]`;
export const DOWNLOAD_BUTTON_ID = "download-button";
// Error messages
export const ERROR_NO_TESTID_OR_LOCATOR =
  "Either testId or locator must be defined";
export const TWO_DECIMAL_NUMBER_REGEX = /^\d+\.?\d{0,2}$/;
