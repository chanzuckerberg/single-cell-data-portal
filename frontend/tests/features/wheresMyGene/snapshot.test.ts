/**
 * Test suite for select filter-related utils.
 */
import { expect, test } from "@playwright/test";
import path from "path";
import {
  checkSourceData,
  deSelectFIlterOption,
  goToWMG,
  selectFIlterOption,
  selectTissueandGeneOption,
} from "../../utils/xgeneUtils";
import { isDevStagingProd } from "tests/utils/helpers";
import { getText } from "tests/utils/selectors";
const CHEVRON_LEFT = '[data-icon="chevron-left"]';

const { describe, skip } = test;

describe("Left side bar", () => {
  //skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test.beforeEach(async ({ page }) => {
    // navigate to url
    await goToWMG(page);
    await selectTissueandGeneOption(page);
  });

  test("Left side bar collapse and expand", async ({ page }) => {
    // click chevron left to collapse the left tab
    await page.locator(CHEVRON_LEFT).click();

    // verify the left tab is collapsed
    expect(await page.getByTestId("add-organism").isVisible()).toBeFalsy();
  });

  test.only("should show a visual change when selecting and de-selecting filter option", async ({
    page,
  }) => {
    await page.pause();
    const elementHandle = await page.locator('[id="heatmap-container-id"]');

    const beforeScreenshotPath = path.resolve("../before-filter.png");
    const afterScreenshotPath = path.resolve("../after-filter.png");
    await elementHandle.screenshot({ path: beforeScreenshotPath });
    await selectFIlterOption(page, "dataset-filter");

    //wait for filter to reflect on the UI
    await page.locator(getText("Loading")).waitFor({ state: "hidden" });
    await elementHandle.screenshot({ path: afterScreenshotPath });

    expect(
      await elementHandle.screenshot({ path: afterScreenshotPath })
    ).toMatchSnapshot(beforeScreenshotPath, { maxDiffPixels: 5 });

    await deSelectFIlterOption(page, "dataset-filter");
  });
});
