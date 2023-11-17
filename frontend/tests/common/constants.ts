import { CELL_TYPE_NAME_LABEL_CLASS_NAME } from "src/views/WheresMyGeneV2/components/HeatMap/components/YAxisChart/constants";

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

export const SKIP_LOGIN =
  process.env.SKIP_LOGIN === "true" || process.env.LOGIN === "false" || false;

export const LOGIN_STATE_FILENAME = "loginState.json";

export const SHARED_LINK = `${TEST_URL}/gene-expression?compare=disease&diseases=MONDO%3A0100096&ethnicities=unknown&sexes=PATO%3A0000383%2CPATO%3A0000384&tissues=UBERON%3A0000178%2CUBERON%3A0002048&genes=DPM1%2CTNMD%2CTSPAN6%2CMALAT1&ver=2`;

export const SHARED_LINK_WITHOUT_COMPARE = `${TEST_URL}/gene-expression?diseases=MONDO%3A0100096&ethnicities=unknown&sexes=PATO%3A0000383%2CPATO%3A0000384&tissues=UBERON%3A0000178%2CUBERON%3A0002048&genes=DPM1%2CTNMD%2CTSPAN6%2CMALAT1&ver=2`;

export const SIMPLE_SHARED_LINK = `${TEST_URL}/gene-expression?tissues=UBERON%3A0000178%2CUBERON%3A0002048&genes=DPM1%2CMALAT1%2CTNMD%2CTSPAN6&ver=2`;

export const SHARED_LINK_NO_FILTER = `${TEST_URL}/gene-expression?tissues=UBERON%3A0000178%2CUBERON%3A0002048&genes=DPM1%2CTNMD%2CTSPAN6&ver=2`;

export const SHARED_LINK_FILTER = `${TEST_URL}/gene-expression?compare=disease&sexes=PATO%3A0000383&tissues=UBERON%3A0000178%2CUBERON%3A0002048&genes=DPM1%2CTNMD%2CTSPAN6&ver=2`;
export const SHARED_LINK_NO_GROUP = `${TEST_URL}/gene-expression?compare=disease&diseases=PATO%3A0000461&tissues=UBERON%3A0000178%2CUBERON%3A0002048&genes=DPM1%2CTNMD%2CTSPAN6&ver=2`;

export const ADD_GENE_BTN = "add-gene-btn";

export const CELL_TYPE_LABELS_ID = CELL_TYPE_NAME_LABEL_CLASS_NAME;
export const TWO_DECIMAL_NUMBER_REGEX = /^\d+\.?\d{0,2}$/;
export const ERROR_NO_TESTID_OR_LOCATOR =
  "Either testId or locator must be defined";
export const GENE_LABELS_ID = "[data-testid^=gene-label-]";

export const downloadPath = "./tests/downloads";
export const HOMO_SAPIENS_TERM_ID = "NCBITaxon:9606";
export const ADD_TISSUE_ID = "tissue-filter";
export const ADD_GENE_ID = "add-gene-btn";
export const GENE_DELETE_BUTTON = "gene-delete-button";
export const SOURCE_DATA_BUTTON_ID = "source-data-button";
export const SOURCE_DATA_LIST_SELECTOR = `[data-testid="source-data-list"]`;
export const DOWNLOAD_BUTTON_ID = "download-button";
export const COMPARE_DROPDOWN_ID = "compare-dropdown";
export const SEX_FILTER_TEST_ID = "sex-filter";
export const PUBLICATION_FILTER_TEST_ID = "publication-filter";
export const DISEASE_FILTER_TEST_ID = "disease-filter";
export const CELL_TYPE_FILTER_TEST_ID = "celltype-filter";
export const SELF_REPORTED_ETHNICITY_FILTER_TEST_ID =
  "self-reported-ethnicity-filter";
