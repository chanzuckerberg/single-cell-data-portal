/**
 * Test suite for select filter-related utils.
 */
import { expect, test } from "@playwright/test";
import { goToWMG, selectAndDeselectOption } from "tests/utils/xgeneUtils";
import { isDevStagingProd } from "tests/utils/helpers";
const CHEVRON_LEFT = '[data-icon="chevron-left"]';

const { describe, skip } = test;

describe("Left side bar", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
  test.beforeEach(async ({ page }) => {
    // navigate to url
    await goToWMG(page);
  });

  test("Left side bar collapse and expand", async ({ page }) => {
    // click chevron left to collapse the left tab
    await page.locator(CHEVRON_LEFT).click();

    // verify the left tab is collapsed
    expect(await page.getByTestId("dataset-filter").isVisible()).toBeFalsy();
  });

  test("Should be able select and de-select options for datasetfilter", async ({
    page,
  }) => {
    await selectAndDeselectOption(page, "dataset-filter");
  });

  test("Should be able select and de-select options for disease filter", async ({
    page,
  }) => {
    await selectAndDeselectOption(page, "disease-filter");
  });

  test("Should be able select and de-select options for Self-Reported Ethnicity filter", async ({
    page,
  }) => {
    await selectAndDeselectOption(page, "self-reported-ethnicity-filter");
  });

  test("Should be able select and de-select options for Self-Reported sex-filter", async ({
    page,
  }) => {
    await selectAndDeselectOption(page, "sex-filter");
  });
});
