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
const geneCanvasId = '[data-zr-dom-id*="zr"]';
describe("SVG download tests", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test(`Should verify SVG download without grouping`, async ({ page }) => {
    // set app state
    await goToWMG(page, SIMPLE_SHARED_LINK);

    const tissues = ["blood"]; // todo: handle multiple tissues
    const fileTypes = ["svg"];
    const folder = subDirectory();
    //download and verify svg file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);

    // verify SVG

    for (let i = 0; i < tissues.length; i++) {
      const cellSnapshot = `${downLoadPath}/${tissues[i]}.png`;
      const geneSnapshot = `${downLoadPath}/gene_${i}.png`;
      // Get the current viewport size
      const currentViewportSize = page.viewportSize();
      // Find the y-axis element
      const yAxisElement = page.locator(getById(`${tissues[i]}-y-axis`));

      // Get the bounding box of the element
      const box = await yAxisElement.boundingBox();
      if (!box) {
        console.error("Element not found");
        return;
      }
      // Set the viewport size to the size of the element
      // some of the elements were missing in the screenshot
      await page.setViewportSize({
        width: box.width * 10,
        height: box.height,
      });

      // Take a screenshot of the y-axis element
      await yAxisElement.screenshot({
        path: cellSnapshot,
      });

      await page
        .locator(geneCanvasId)
        .nth(0)
        .screenshot({ path: geneSnapshot });
      // Revert the viewport size to its original value
      if (currentViewportSize) {
        await page.setViewportSize(currentViewportSize);
      }

      //I need the page to revert to original viewport before  proceeding
      setTimeout(async () => {
        await compareSvg(
          page,
          geneSnapshot,
          geneSnapshot,
          `${folder}/${tissues[i]}.svg`,
          folder
        );
      }, 3000);
    }
  });

  test(`Should verify SVG download with grouping`, async ({ page }) => {
    // set app state
    await goToWMG(page, SHARED_LINK);

    const tissues = ["blood"]; // todo: handle multiple tissues
    const fileTypes = ["svg"];
    const folder = subDirectory();
    //download and verify svg file
    await downloadAndVerifyFiles(page, fileTypes, tissues, folder);

    // verify SVG

    for (let i = 0; i < tissues.length; i++) {
      const cellSnapshot = `${downLoadPath}/${tissues[i]}.png`;
      const geneSnapshot = `${downLoadPath}/gene_${i}.png`;
      // Get the current viewport size
      const currentViewportSize = page.viewportSize();
      // Find the y-axis element
      const yAxisElement = page.locator(getById(`${tissues[i]}-y-axis`));

      // Get the bounding box of the element
      const box = await yAxisElement.boundingBox();
      if (!box) {
        console.error("Element not found");
        return;
      }
      // Set the viewport size to the size of the element
      // some of the elements were missing in the screenshot
      await page.setViewportSize({
        width: box.width * 10,
        height: box.height,
      });

      // Take a screenshot of the y-axis element
      await yAxisElement.screenshot({
        path: cellSnapshot,
      });

      await page
        .locator(geneCanvasId)
        .nth(0)
        .screenshot({ path: geneSnapshot });
      // Revert the viewport size to its original value
      if (currentViewportSize) {
        await page.setViewportSize(currentViewportSize);
      }

      //I need the page to revert to original viewport before  proceeding
      setTimeout(async () => {
        await compareSvg(
          page,
          geneSnapshot,
          geneSnapshot,
          `${folder}/${tissues[i]}.svg`,
          folder
        );
      }, 3000);
    }
  });
  test.afterAll(async () => {
    deleteDownloadedFiles(`./tests/downloads`);
  });
});
