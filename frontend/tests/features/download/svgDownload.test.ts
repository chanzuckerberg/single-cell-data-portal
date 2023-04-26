import { expect, test } from "@playwright/test";
import {
  deleteCsvDownloads,
  downloadCsv,
  downloadSvg,
  getCsvMetadata,
  goToWMG,
  selectTissueAndGeneOption,
} from "../../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";
import fs from "fs";
const BLOOD_TISSUE_COUNT =
  '[data-testid="cell-type-labels-blood"] [data-testid="cell-type-label-count"]';
const { describe, skip } = test;
describe("Csv download", () => {
  // skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test.beforeEach(async ({ page }) => {
    // navigate to gene expression page
    await goToWMG(page);
    //select tissue and gene
    await selectTissueAndGeneOption(page);
  });

  test(`Verify metadata and header displayed on csv file with  applied`, async ({
    page,
  }) => {
    // // to differentiate file in each test run
    const randomNumber: number = Math.floor(Math.random() * 90000) + 10000;
    const fileFactor: string = randomNumber.toString();
    const tissue = "blood";
    //download and extract the csv file
    await downloadSvg(page, fileFactor);
    const svgFile = fs.readFileSync(
      `./tests/download/${fileFactor}/${tissue}.svg`,
      "utf-8"
    );
    //  await page.pause();
    await page.setContent(svgFile);

    // Use page.evaluate() to interact with the SVG elements
    const viewBox = await page.evaluate(() => {
      const svg = document.querySelector("svg");

      // Count the number of elements with id "cell-name-label-group"
      if (svg) {
        console.log("the endddd.d.d..d.d");
        return svg.querySelectorAll("#cell-name-label-group").length;
      }
      return 0;
    });

    console.log(viewBox);
  });
  test.only(`Verify metadata and header displayed on csv file with applied second`, async ({
    page,
  }) => {
    const randomNumber: number = Math.floor(Math.random() * 90000) + 10000;
    const fileFactor: string = randomNumber.toString();
    const tissue = "blood";

    // Download and extract the csv file
    await downloadSvg(page, fileFactor);

    const svgFile = fs.readFileSync(
      `./tests/download/${fileFactor}/${tissue}.svg`,
      "utf-8"
    );
    await page.setContent(svgFile);

    // Use page.evaluate() to interact with the SVG elements
    const labels = await page.evaluate(() => {
      const svg = document.querySelector("svg");

      // Extract the text content of elements with id "cell-name-label-group"
      if (svg) {
        const labelElements = Array.from(
          svg.querySelectorAll("#cell-name-label-group")
        );
        return labelElements.map(
          (el) => el.textContent && el.textContent.trim()
        );
      }
      return [];
    });

    console.log(labels);
  });

  test.afterAll(async () => {
    //delete csv
    // deleteCsvDownloads(`./tests/download`);
  });
});
