/**
 * Test suite for select filter-related utils.
 */
import { test } from "@playwright/test";
import {
  goToWMG,
  selectFilterOption,
  selectTissueandGeneOption,
} from "tests/common/utils";
import { isDevStagingProd } from "tests/utils/helpers";

const { describe, skip } = test;

describe("Left side bar", () => {
  //skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test.beforeEach(async ({ page }) => {
    // navigate to url
    await page.pause();
    await goToWMG(page);
    await selectTissueandGeneOption(page);
  });

  test.only("Should check source after a filter is slected ", async ({
    page,
  }) => {
    await selectFilterOption(page, "sex-filter");
  });
});
