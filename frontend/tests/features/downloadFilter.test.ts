import { expect, test } from "@playwright/test";
import {
  downloadCsv,
  getCsvHeaders,
  getCsvMetadata,
  getFilterApplied,
  goToWMG,
  selectFilterOption,
  selectTissueAndGeneOption,
} from "../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";

const EXPECTED_HEADER = [
  "Tissue",
  "Cell Type",
  "Cell Count",
  "Tissue Composition",
  "Gene Symbol",
  "Expression",
  '"Expression',
  ' Scaled"',
  "Number of Cells Expressing Genes",
];
const EXPECTED_METADATA = [
  "# We regularly expand our single cell data corpus to improve results. Downloaded data and figures may differ in the future.",
  "# Dataset Filter Values: No selection",
  "# Disease Filter Values: No selection",
  "# Self-Reported Ethnicity Filter Values: No selection",
  "# Sex Filter Values: No selection",
  "# Organism Filter Value: Homo sapiens",
];
const { describe, skip } = test;
describe("Csv download", () => {
  // skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test.beforeEach(async ({ page }) => {
    // navigate to gene expression page
    await goToWMG(page);

    //select tissue and gene
    await selectTissueAndGeneOption(page);
  });

  test.only("Verify metadata displayed on csv file with filter applied", async ({
    page,
  }) => {
    //select filter
    await selectFilterOption(page, "dataset-filter");
    const nur = await getFilterApplied(page, "dataset-filter");

    //download and extract the csv file
    // await downloadCsv(page);
  });
  test("Verify headers displayed on csv file", async () => {
    // put all the headers in an array
    const HEADERS = await getCsvHeaders("blood");

    //verify all the headers are present in the csv
    expect(HEADERS[0]).toEqual(expect.arrayContaining(EXPECTED_HEADER));
  });
});
