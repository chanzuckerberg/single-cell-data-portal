import { expect, test } from "@playwright/test";
import {
  downloadCsv,
  getCsvHeaders,
  getCsvMetadata,
  goToWMG,
  selectTissueAndGeneOption,
} from "../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";

const { describe, skip } = test;
describe("Csv download", () => {
  // skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test.beforeEach(async ({ page }) => {
    // navigate to gene expression page
    await goToWMG(page);

    //select tissue and gene
    await selectTissueAndGeneOption(page);
  });

  test("Verify metadata displayed on csv file", async ({ page }) => {
    // await page.pause();
    await downloadCsv(page);
    const META_DATA = await getCsvMetadata("blood");

    // verify the date is valid
    const dateString = META_DATA[0].substring(14);
    const date = new Date(dateString);
    // Check if the resulting date is valid
    expect(!isNaN(date.getTime())).toBe(true);
    //expect ( META_DATA[2]).toMatchText("https://localhost:3000/gene-expression")
    const EXPECTED_METADATA = [
      "# We regularly expand our single cell data corpus to improve results. Downloaded data and figures may differ in the future.",
      "# Dataset Filter Values: No selection",
      "# Disease Filter Values: No selection",
      "# Self-Reported Ethnicity Filter Values: No selection",
      "# Sex Filter Values: No selection",
      "# Organism Filter Value: Homo sapiens",
    ];

    expect(META_DATA).toEqual(expect.arrayContaining(EXPECTED_METADATA));
  });
  test.only("Verify headers displayed on csv file", async ({ page }) => {
    // await page.pause();
    await downloadCsv(page);
    const HEADERS = await getCsvHeaders("blood");
    const expectedHeaders = [
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

    expect(HEADERS[0]).toEqual(expect.arrayContaining(expectedHeaders));
  });
});
