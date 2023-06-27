import { expect, test } from "@playwright/test";
import {
  checkPlotSize,
  checkSourceData,
  conditionallyRunTests,
  deSelectSecondaryFilterOption,
  goToWMG,
  selectSecondaryFilterOption,
  selectTissueAndGeneOption,
} from "../../utils/wmgUtils";
import { tryUntil } from "tests/utils/helpers";
import {
  COLOR_SCALE_TOOLTIP_TEXT,
  GROUP_BY_TOOLTIP_TEXT,
  SORT_CELL_TYPES_TOOLTIP_TEXT,
  SORT_GENES_TOOLTIP_TEXT,
} from "src/views/WheresMyGene/common/constants";

const SIDE_BAR_TOGGLE_BUTTON_ID = "side-bar-toggle-button";

const { describe } = test;

describe("Left side bar", () => {
  conditionallyRunTests();

  test("Left side bar collapse and expand", async ({ page }) => {
    // navigate to gene expression page
    await goToWMG(page);

    //select tissue and gene
    await selectTissueAndGeneOption(page);

    // click chevron left to collapse the left tab
    await page.getByTestId(SIDE_BAR_TOGGLE_BUTTON_ID).click();

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
          await selectSecondaryFilterOption(page, filterOption);

          // check the count of source data displayed after adding a filter
          const countAfterFilter = await checkSourceData(page);

          //check plot height after adding a filter
          const plotSizeAfterFilter = await checkPlotSize(page);

          // verify source data loading some data
          expect(countBeforeFilter).toBeGreaterThan(0);

          // verify source data changed after filter is applied
          expect(countBeforeFilter === countAfterFilter).toBeFalsy();

          //verify data plot data loading some data
          expect(plotSizeBeforeFilter).toBeGreaterThan(0);

          //verify data plot data changed after filter was applied
          expect(plotSizeBeforeFilter === plotSizeAfterFilter).toBeFalsy();

          //uncheck filter
          await deSelectSecondaryFilterOption(page, filterOption);
        },
        { page }
      );
    });
  });

  test("Left side bar tooltips", async ({ page }) => {
    // navigate to gene expression page
    await goToWMG(page);

    //select tissue and gene
    await selectTissueAndGeneOption(page);

    // Group By tooltip
    await page.getByTestId("group-by-tooltip-icon").hover();
    expect(page.getByText(GROUP_BY_TOOLTIP_TEXT)).toBeTruthy();

    // Color Scale tooltip
    await page.getByTestId("color-scale-tooltip-icon").hover();
    expect(page.getByText(COLOR_SCALE_TOOLTIP_TEXT)).toBeTruthy();

    // Sort Cell Type tooltip
    await page.getByTestId("sort-cell-types-tooltip-icon").hover();
    expect(page.getByText(SORT_CELL_TYPES_TOOLTIP_TEXT)).toBeTruthy();

    // Sort Gene tooltip
    await page.getByTestId("sort-genes-tooltip-icon").hover();
    expect(page.getByText(SORT_GENES_TOOLTIP_TEXT)).toBeTruthy();
  });
});
