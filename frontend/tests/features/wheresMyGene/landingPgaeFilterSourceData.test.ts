/**
 * Test suite for select filter-related utils.
 */
import { test } from "@playwright/test";
import {
  goToWMG,
  selectFilterOption,
  selectTissueandGeneOption,
} from "tests/utils/xgeneUtils";
import { isDevStagingProd } from "tests/utils/helpers";

const { describe, skip } = test;

describe("Left side bar", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test("Should check source after a filter is selected ", async ({ page }) => {
    await goToWMG(page);
    await selectTissueandGeneOption(page);
    await selectFilterOption(page, "sex-filter");
  });
});
