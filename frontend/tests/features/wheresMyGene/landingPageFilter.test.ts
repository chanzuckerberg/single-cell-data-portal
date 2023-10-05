import { expect, test, Page } from "@playwright/test";
import {
  WMG_WITH_SEEDED_GENES,
  checkPlotSize,
  checkSourceData,
  conditionallyRunTests,
  deSelectSecondaryFilterOption,
  goToWMG,
  selectSecondaryFilterOption,
  waitForHeatmapToRender,
} from "tests/utils/wmgUtils";
import { goToPage, tryUntil } from "tests/utils/helpers";
import {
  COLOR_SCALE_TOOLTIP_TEXT,
  GROUP_BY_TOOLTIP_TEXT,
  SORT_CELL_TYPES_TOOLTIP_TEXT,
  SORT_GENES_TOOLTIP_TEXT,
} from "src/views/WheresMyGene/common/constants";

const SIDE_BAR_TOGGLE_BUTTON_ID = "side-bar-toggle-button";
const CELL_TYPE_FILTER = "naive B cell";

const { describe } = test;

describe("Left side bar", () => {
  /**
   * TODO(thuang): Remove forceRun when all WMG e2e tests are enabled.
   * `forceRun` is just to incrementally add tests back in the meantime
   */
  conditionallyRunTests({ forceRun: true });

  test("Left side bar collapse and expand", async ({ page }) => {
    await goToWMG(page);

    // click chevron left to collapse the left tab
    await page.getByTestId(SIDE_BAR_TOGGLE_BUTTON_ID).click();

    // verify the left tab is collapsed
    expect(await page.getByTestId("add-organism").isVisible()).toBeFalsy();
  });

  [
    "dataset-filter",
    "disease-filter",
    "self-reported-ethnicity-filter",
    "sex-filter",
  ].forEach((filterOption) => {
    test(`Should be able select and de-select options for ${filterOption} filter`, async ({
      page,
    }) => {
      await goToPage(WMG_WITH_SEEDED_GENES.URL, page);

      await waitForHeatmapToRender(page);
      await tryUntil(
        async () => {
          // check the count of source data displayed before adding a filter
          const countBeforeFilter = await checkSourceData(page);

          // verify source data loading some data
          expect(countBeforeFilter).toBeGreaterThan(0);

          // check plot height before adding a filter
          const plotSizeBeforeFilter = await checkPlotSize(page);

          // verify data plot data loading some data
          expect(plotSizeBeforeFilter).toBeGreaterThan(0);

          // select a filter
          await selectSecondaryFilterOption(page, filterOption);

          await tryUntil(
            async () => {
              // check the count of source data displayed after adding a filter
              const countAfterFilter = await checkSourceData(page);

              //check plot height after adding a filter
              const plotSizeAfterFilter = await checkPlotSize(page);

              // verify source data changed after filter is applied
              expect(countBeforeFilter === countAfterFilter).toBeFalsy();

              // verify data plot data changed after filter was applied
              expect(plotSizeBeforeFilter === plotSizeAfterFilter).toBeFalsy();
            },
            { page }
          );

          // uncheck filter
          await deSelectSecondaryFilterOption(page, filterOption);
        },
        {
          page,
          /**
           * (thuang): Give up after N times, because the app state might not
           * be recoverable at this point
           */
          maxRetry: 3,
        }
      );
    });
  });

  test("Left side bar tooltips", async ({ page }) => {
    await goToPage(WMG_WITH_SEEDED_GENES.URL, page);

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

  [
    "dataset-filter",
    "disease-filter",
    "publication-filter",
    "tissue-filter",
  ].forEach((testId) => {
    test(`Ensure that cell type filter cross-filters with ${testId}`, async ({
      page,
    }) => {
      await goToPage(WMG_WITH_SEEDED_GENES.URL, page);

      await waitForHeatmapToRender(page);

      const numberOfRecordsBefore = await countRecords(page, testId);

      await page.getByRole("combobox").first().click();
      await page.getByRole("option", { name: CELL_TYPE_FILTER }).click();
      await page.keyboard.press("Escape");

      const numberOfRecordsAfter = await countRecords(page, testId);

      expect(numberOfRecordsBefore).toBeGreaterThan(numberOfRecordsAfter);
    });
  });

  async function countRecords(page: Page, testId: string) {
    const optionLocator = page.getByRole("option");

    let count = 0;

    await tryUntil(
      async () => {
        await page.getByTestId(testId).getByRole("button").click();

        count = await optionLocator.count();

        await page.keyboard.press("Escape");

        expect(count).toBeGreaterThan(0);
      },
      { page }
    );

    return count;
  }
});
