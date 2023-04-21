import { expect, test } from "@playwright/test";
import {
  checkPlotSize,
  checkSourceData,
  deSelectFilterOption,
  goToWMG,
  selectFilterOption,
  selectTissueAndGeneOption,
} from "../../utils/wmgUtils";
import { isDevStagingProd, tryUntil } from "tests/utils/helpers";
const CHEVRON_LEFT = '[data-icon="chevron-left"]';

const { describe, skip } = test;

describe("Left side bar", () => {
  //skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test("Left side bar collapse and expand", async ({ page }) => {
    // navigate to gene expression page
    await goToWMG(page);

    //select tissue and gene
    await selectTissueAndGeneOption(page);

    // click chevron left to collapse the left tab
    await page.locator(CHEVRON_LEFT).click();

    // verify the left tab is collapsed
    expect(await page.getByTestId("add-organism").isVisible()).toBeFalsy();
  });

  [
    ["dataset-filter"],
    ["disease-filter"],
    ["self-reported-ethnicity-filter"],
    ["sex-filter"],
  ].forEach(([filterOption]) => {
    test(`Should be able select and de-select options for ${filterOption} filter`, async ({
      page,
    }) => {
      // navigate to gene expression page
      await goToWMG(page);

      //select tissue and gene
      await selectTissueAndGeneOption(page);

      await tryUntil(
        async () => {
          // check the count of source data displayed before adding a filter
          const countBeforeFilter = await checkSourceData(page);

          //check plot height before adding a filter
          const plotSizeBeforeFilter = await checkPlotSize(page);

          //select a filter
          await selectFilterOption(page, filterOption);

          // check the count of source data displayed after adding a filter
          const countAfterFilter = await checkSourceData(page);

          //check plot height after adding a filter
          const plotSizeAfterFilter = await checkPlotSize(page);

          // verify source  data loading some data
          expect(countBeforeFilter).toBeGreaterThan(0);
          // verify source  data changed after filter is applied
          expect(countBeforeFilter === countAfterFilter).toBeFalsy();

          //verify data plot data loading some data
          expect(plotSizeBeforeFilter).toBeGreaterThan(0);

          //verify data plot data changed after filter was  applied
          expect(plotSizeBeforeFilter === plotSizeAfterFilter).toBeFalsy();
          //uncheck filter
          await deSelectFilterOption(page, filterOption);
        },
        { page }
      );
    });
  });
});
