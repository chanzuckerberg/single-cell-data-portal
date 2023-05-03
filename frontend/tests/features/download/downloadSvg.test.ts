import { test } from "@playwright/test";
import { goToWMG } from "../../utils/wmgUtils";
import {
  SHARED_LINK,
  SIMPLE_SHARED_LINK,
  downLoadPath,
} from "tests/common/constants";
import {
  compareSvg,
  deleteDownloadedFiles,
  downloadAndVerifyFiles,
  subDirectory,
} from "tests/utils/downloadUtils";
import { isDevStagingProd } from "tests/utils/helpers";
import { getById } from "tests/utils/selectors";

const { describe, skip } = test;
const geneCanvasId = '[data-zr-dom-id="zr_100000"]';
describe("SVG download tests", () => {
  //skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test.only(`Should verify SVG download without grouping`, async ({ page }) => {
    // set app state
    await goToWMG(page, SIMPLE_SHARED_LINK);

    const tissues = ["blood", "lung"];
    const fileTypes = ["svg"];
    const folder = subDirectory();
    //download and verify svg file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);

    // verify SVG

    for (let i = 0; i < tissues.length; i++) {
      const cellSnapshot = `${downLoadPath}/${tissues[i]}.png`;
      const geneSnapshot = `${downLoadPath}/gene_${i}.png`;
      await page
        .locator(getById(`${tissues[i]}-y-axis`))
        .screenshot({ path: cellSnapshot });
      await page
        .locator(geneCanvasId)
        .nth(0)
        .screenshot({ path: geneSnapshot });
      await compareSvg(geneSnapshot, geneSnapshot, `${tissues[i]}.svg`);
    }
  });

  test(`Should verify SVG download with grouping`, async ({ page }) => {
    // set app state
    await goToWMG(page, SIMPLE_SHARED_LINK);

    const tissues = ["blood", "lung"];
    const fileTypes = ["svg"];
    const folder = subDirectory();
    //download and verify svg file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);

    for (let i = 0; i < tissues.length; i++) {
      const cellSnapshot = `${downLoadPath}/${tissues[i]}.png`;
      const geneSnapshot = `${downLoadPath}/gene_${i}.png`;
      await page
        .locator(getById(`${tissues[i]}-y-axis`))
        .screenshot({ path: cellSnapshot });
      await page
        .locator(geneCanvasId)
        .nth(0)
        .screenshot({ path: geneSnapshot });
      await compareSvg(geneSnapshot, geneSnapshot, `${tissues[i]}.svg`);
    }
  });
  // test.afterAll(async () => {
  //   //delete csv
  //   deleteDownloadedFiles(`./tests/downloads`);
  // });
});
