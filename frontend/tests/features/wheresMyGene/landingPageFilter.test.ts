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
  selectTissueAndGeneOption,
} from "../../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";
const CHEVRON_LEFT = '[data-icon="chevron-left"]';

const { describe, skip } = test;

describe("Left side bar", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test("Left side bar collapse and expand", async ({ page }) => {
    await goToWMG(page);
    await selectTissueAndGeneOption(page);
    // click chevron left to collapse the left tab
    await page.locator(CHEVRON_LEFT).click();

    // verify the left tab is collapsed
    expect(await page.getByTestId("add-organism").isVisible()).toBeFalsy();
  });

  test.only("Should be able select and de-select options for datasetfilter", async ({
    page,
  }) => {
    await goToWMG(page);
    await selectTissueAndGeneOption(page);
    const countBeforeFIlter = await checkSourceData(page);
    const plotSizeBeforeFIlter = await checkPlotSize(page);
    await selectFIlterOption(page, "dataset-filter");
    const countAfterFIlter = await checkSourceData(page);
    const plotSizeAfterFIlter = await checkPlotSize(page);
    expect(countBeforeFIlter).toBeGreaterThan(0);
    expect(countBeforeFIlter === countAfterFIlter).toBeFalsy();
    expect(plotSizeBeforeFIlter).toBeGreaterThan(0);
    expect(plotSizeBeforeFIlter === plotSizeAfterFIlter).toBeFalsy();
    await deSelectFIlterOption(page, "dataset-filter");
  });

  test("Should be able select and de-select options for disease filter", async ({
    page,
  }) => {
    await goToWMG(page);
    await selectTissueAndGeneOption(page);
    const countBeforeFIlter = await checkSourceData(page);
    const plotSizeBeforeFIlter = await checkPlotSize(page);
    await selectFIlterOption(page, "disease-filter");
    const countAfterFIlter = await checkSourceData(page);
    const plotSizeAfterFIlter = await checkPlotSize(page);
    
    expect(countBeforeFIlter).toBeGreaterThan(0);
    expect(countBeforeFIlter === countAfterFIlter).toBeFalsy();
    expect(plotSizeBeforeFIlter).toBeGreaterThan(0);
    expect(plotSizeBeforeFIlter === plotSizeAfterFIlter).toBeFalsy();
    await deSelectFIlterOption(page, "disease-filter");
  });

  test("Should be able select and de-select options for Self-Reported Ethnicity filter", async ({
    page,
  }) => {
    await goToWMG(page);
    await selectTissueAndGeneOption(page);
    const countBeforeFIlter = await checkSourceData(page);
    const plotSizeBeforeFIlter = await checkPlotSize(page);
    await selectFIlterOption(page, "self-reported-ethnicity-filter");
    const countAfterFIlter = await checkSourceData(page);
    const plotSizeAfterFIlter = await checkPlotSize(page);
    expect(countBeforeFIlter).toBeGreaterThan(0);
    expect(countBeforeFIlter === countAfterFIlter).toBeFalsy();
    expect(plotSizeBeforeFIlter).toBeGreaterThan(0);
    expect(plotSizeBeforeFIlter === plotSizeAfterFIlter).toBeFalsy();
    await deSelectFIlterOption(page, "self-reported-ethnicity-filter");
  });

  test("Should be able select and de-select options for Self-Reported sex-filter", async ({
    page,
  }) => {
    await goToWMG(page);
    await selectTissueAndGeneOption(page);

    const countBeforeFIlter = await checkSourceData(page);
    const plotSizeBeforeFIlter = await checkPlotSize(page);
    await selectFIlterOption(page, "sex-filter");
    const countAfterFIlter = await checkSourceData(page);
    const plotSizeAfterFIlter = await checkPlotSize(page);
    expect(countBeforeFIlter).toBeGreaterThan(0);
    expect(countBeforeFIlter === countAfterFIlter).toBeFalsy();
    expect(plotSizeBeforeFIlter).toBeGreaterThan(0);
    expect(plotSizeBeforeFIlter === plotSizeAfterFIlter).toBeFalsy();
    await deSelectFIlterOption(page, "sex-filter");
  });
});
