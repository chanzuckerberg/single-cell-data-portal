import { expect, test } from "@playwright/test";
import {
  deleteCsvFile,
  downloadCsv,
  getCsvHeaders,
  getCsvMetadata,
  goToWMG,
  selectFilterOption,
  selectTissueAndGeneOption,
  verifyMetadata,
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

const { describe, skip } = test;
describe("Csv download", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test.beforeEach(async ({ page }) => {
    // navigate to gene expression page
    await goToWMG(page);
    //select tissue and gene
    await selectTissueAndGeneOption(page);
  });
  [
    ["dataset-filter"],
    ["disease-filter"],
    ["self-reported-ethnicity-filter"],
    ["sex-filter"],
  ].forEach(([filter]) => {
    test(`Verify metadata and header displayed on csv file with ${filter} applied`, async ({
      page,
      context,
    }) => {
      // to differentiate file in each test run
      const randomNumber: number = Math.floor(Math.random() * 90000) + 10000;
      const fileFactor: string = randomNumber.toString();
      //select filter
      await selectFilterOption(page, filter);
      //download and extract the csv file
      await downloadCsv(page, fileFactor);

      //put all the meta data in an array
      const data = await getCsvMetadata("blood", fileFactor);

      // put all the headers in an array
      const headers = await getCsvHeaders("blood", fileFactor);

      //verify meta data
      await verifyMetadata(page, filter, data, context);

      //verify all the headers are present in the csv
      expect(headers[0]).toEqual(expect.arrayContaining(EXPECTED_HEADER));
    });
  });
  test("Verify metadata displayed on csv file with no filter", async ({
    page,
    context,
  }) => {
    // to differentiate file in each test run
    const randomNumber: number = Math.floor(Math.random() * 90000) + 10000;
    const fileFactor: string = randomNumber.toString();
    //download and extract the csv file
    await downloadCsv(page, fileFactor);

    // put all the headers in an array
    const headers = await getCsvHeaders("blood", fileFactor);

    //put all the meta data in an array
    const data = await getCsvMetadata("blood", fileFactor);

    //verify meta data
    await verifyMetadata(page, "no-filter", data, context);

    //verify all the headers are present in the csv
    expect(headers[0]).toEqual(expect.arrayContaining(EXPECTED_HEADER));
  });
  test.afterAll(async () => {
    //delete csv
    deleteCsvFile(`./tests/utils/download`);
  });
});
