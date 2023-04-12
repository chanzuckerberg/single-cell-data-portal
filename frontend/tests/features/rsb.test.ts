/**
 * Test suite for select filter-related utils.
 */
import { expect, test } from "@playwright/test";
import {
  checkPlotSize,
  checkSourceData,
  deSelectFIlterOption,
  goToWMG,
  selectFIlterOption,
  selectTissueandGeneOption,
} from "../../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";
const CHEVRON_LEFT = '[data-icon="chevron-left"]';

const { describe, skip } = test;

describe("Right side bar", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test.beforeEach(async ({ page }) => {
    // navigate to url
    await goToWMG(page);
  });

  test("Left side bar collapse and expand", async ({ page }) => {
    // click chevron left to collapse the left tab
    await page.locator(CHEVRON_LEFT).click();

    // verify the left tab is collapsed
    expect(await page.getByTestId("add-organism").isVisible()).toBeFalsy();
  });
});
