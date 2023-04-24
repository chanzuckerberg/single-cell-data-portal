import { expect, test } from "@playwright/test";
import {
  deleteCsvFile,
  downloadCsv,
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
const BLOOD_TISSUE_COUNT =
  '[data-testid="cell-type-labels-blood"] [data-testid="cell-type-label-count"]';
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
    ["no-filter"],
  ].forEach(([filter]) => {
    test.only(`Verify metadata and header displayed on csv file with ${filter} applied`, async ({
      page,
    }) => {
      // to differentiate file in each test run
      const randomNumber: number = Math.floor(Math.random() * 90000) + 10000;
      const fileFactor: string = randomNumber.toString();

      if (filter !== "no-filter") {
        //select filter
        await selectFilterOption(page, filter);
      }

      //download and extract the csv file
      await downloadCsv(page, fileFactor);
      const metadata = await getCsvMetadata("blood", fileFactor);

      // extract the headers and data arrays from the metadata object
      // put all the headers in an array
      const headers = metadata.headers;
      const data = metadata.data;

      //get number of element in csv
      const csvElementsCount = metadata.rowCount;

      //get number of element displayed in ui
      const uiElementsCount = await page.locator(BLOOD_TISSUE_COUNT).count();

      //verify the number of element in teh csv
      expect(csvElementsCount).toEqual(uiElementsCount * 3);

      //verify meta data
      await verifyMetadata(page, filter, headers);

      //verify all the headers are present in the csv
      expect(data[0]).toEqual(expect.arrayContaining(EXPECTED_HEADER));
    });
  });
  test.afterAll(async () => {
    //delete csv
    deleteCsvFile(`./tests/utils/download`);
  });
});