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
    const elementHandle = await page.locator('[id="heatmap-container-id"]');

    const beforeScreenshotPath = path.resolve("../before-screenshot.png");
    const afterScreenshotPath = path.resolve("../after-screenshot.png");

    await elementHandle.screenshot({ path: beforeScreenshotPath });
    await selectFIlterOption(page, "dataset-filter");
    await elementHandle.screenshot({ path: afterScreenshotPath });

    expect(
      await elementHandle.screenshot({ path: afterScreenshotPath })
    ).toMatchSnapshot(beforeScreenshotPath);

    await deSelectFIlterOption(page, "dataset-filter");
  });
});
