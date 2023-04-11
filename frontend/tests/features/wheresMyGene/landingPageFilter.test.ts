/**
 * Test suite for select filter-related utils.
 */
import { expect, test } from "@playwright/test";
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
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");
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

  test.only("Should be able select and de-select options for datasetfilter", async ({
    page,
  }) => {
    const countBeforeFIlter = await checkSourceData(page);
    await selectFIlterOption(page, "dataset-filter");
    const countAfterFIlter = await checkSourceData(page);
    expect(countBeforeFIlter).toBeGreaterThan(0);
    expect(countBeforeFIlter === countAfterFIlter).toBeFalsy();

    await deSelectFIlterOption(page, "dataset-filter");
  });

  test.only("Should be able select and de-select options for disease filter", async ({
    page,
  }) => {
    const countBeforeFIlter = await checkSourceData(page);
    await selectFIlterOption(page, "disease-filter");
    const countAfterFIlter = await checkSourceData(page);
    expect(countBeforeFIlter).toBeGreaterThan(0);
    expect(countBeforeFIlter === countAfterFIlter).toBeFalsy();
    await deSelectFIlterOption(page, "disease-filter");
  });

  test.only("Should be able select and de-select options for Self-Reported Ethnicity filter", async ({
    page,
  }) => {
    const countBeforeFIlter = await checkSourceData(page);
    await selectFIlterOption(page, "self-reported-ethnicity-filter");
    const countAfterFIlter = await checkSourceData(page);
    expect(countBeforeFIlter).toBeGreaterThan(0);
    expect(countBeforeFIlter === countAfterFIlter).toBeFalsy();
    await deSelectFIlterOption(page, "self-reported-ethnicity-filter");
  });

  test.only("Should be able select and de-select options for Self-Reported sex-filter", async ({
    page,
  }) => {
    const countBeforeFIlter = await checkSourceData(page);
    await selectFIlterOption(page, "sex-filter");
    const countAfterFIlter = await checkSourceData(page);
    expect(countBeforeFIlter).toBeGreaterThan(0);
    expect(countBeforeFIlter === countAfterFIlter).toBeFalsy();
    await deSelectFIlterOption(page, "sex-filter");
  });
});
