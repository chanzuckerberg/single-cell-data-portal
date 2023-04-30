import { test } from "@playwright/test";
import { goToWMG, selectTissueAndGeneOption } from "../../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";
import { downloadPng, deleteCsvDownloads } from "tests/utils/downloadUtils";
import pixelmatch from "pixelmatch";
import fs from "fs";
import { PNG } from "pngjs";

const { describe, skip } = test;
describe("CSV download tests", () => {
  test.beforeEach(async ({ page }) => {
    // navigate to gene expression page
    await goToWMG(page);
    //select tissue and gene
    await selectTissueAndGeneOption(page);
  });
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test.only(`Should verify Png for blood tissue with no filter applied`, async ({
    page,
  }) => {
    // generate sub folder
    const randomNumber: number = Math.floor(Math.random() * 90000) + 10000;
    const subDirectory: string = randomNumber.toString();
    //download and verify csv file
    await downloadPng(page, subDirectory);
    // Capture the actual screenshot and compare it with the expected screenshot
    const expectedPng = PNG.sync.read(
      fs.readFileSync(`./tests/fixtures/download/blood.png`)
    );
    const actualPng = PNG.sync.read(
      fs.readFileSync(`./tests/download/${subDirectory}/blood.png`)
    );
    const { width, height } = expectedPng;
    const diff = new PNG({ width, height });

    const mismatchedPixels = pixelmatch(
      expectedPng.data,
      actualPng.data,
      diff.data,
      width,
      height,
      { threshold: 0.1 }
    );

    // Log the number of mismatched pixels
    console.log(`Number of mismatched pixels: ${mismatchedPixels}`);
  });

  test.afterAll(async () => {
    //delete csv
    deleteCsvDownloads(`./tests/download`);
  });
});
