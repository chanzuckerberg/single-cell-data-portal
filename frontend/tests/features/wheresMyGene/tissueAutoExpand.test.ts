import { expect, test, Page } from "@playwright/test";
import { toInteger } from "lodash";
import { collapseTissue, expandTissue } from "tests/utils/helpers";
import { conditionallyRunTests, goToWMG } from "tests/utils/wmgUtils";

const TISSUE_NODE_TEST_ID = "tissue-name";
const TISSUE_FILTER_TEST_ID = "tissue-filter";
const CELL_TYPE_FILTER = "B cell";
const FILTERED_TISSUES = ["abdomen", "axilla", "brain"];

const { describe } = test;

describe("WMG tissue auto-expand", () => {
  conditionallyRunTests({ forceRun: true });
  /**
   * Click the first tissue to expand, filter two tissues from left
   * panel. Expect only two tissues present and both expanded
   * automatically.
   */
  test("Filter tissue auto expansion", async ({ page }) => {
    await loadPageAndTissues(page);
    await expandTissue(page, FILTERED_TISSUES[0]);
    await filterTissues(page, FILTERED_TISSUES);
    await checkTissues(page, FILTERED_TISSUES, FILTERED_TISSUES);
  });
  /**
   * Filter first two tissues from the left panel, collapse the first
   * tissue, add third tissue from the left panel. Expect the first tissue
   * to stay collapsed, and the other two tissues expanded.
   */
  test("Filter tissue auto expansion - respect user collapsed tissue state", async ({
    page,
  }) => {
    await loadPageAndTissues(page);
    await filterTissues(page, FILTERED_TISSUES.slice(0, 2));
    await collapseTissue(page, FILTERED_TISSUES[0]);
    await filterTissues(page, FILTERED_TISSUES.slice(2, 3));
    await checkTissues(page, FILTERED_TISSUES, FILTERED_TISSUES.slice(1, 3));
  });
  /**
   * Filter two tissues, expect tissues to be expanded. Remove both tissues
   * from filter. Expect showing all tissues in collapsed state.
   */
  test("Filter tissue auto expansion - exit filter tissues mode", async ({
    page,
  }) => {
    await loadPageAndTissues(page);
    await filterTissues(page, FILTERED_TISSUES);
    await checkTissues(page, FILTERED_TISSUES, FILTERED_TISSUES);
    await filterTissues(page, FILTERED_TISSUES);
    await checkTissues(page, [], []);
  });
  /**
   * Filter two tissues, expect tissues to be expanded. Expect tissue list
   * to only show tissues that contain "B cell" and automatically expanded.
   * Check that only 'B Cell' cells are visible under expanded tissues.
   * When cell type filter is removed, expect all cells to be shown under
   * expanded tissues.
   */
  test("Filter tissue auto expansion - filter cell type 'B cell'", async ({
    page,
  }) => {
    await loadPageAndTissues(page);
    await filterTissues(page, FILTERED_TISSUES);
    await filterCellType(page);
    await checkTissues(page, FILTERED_TISSUES, FILTERED_TISSUES);
    await checkCellTypes(page);
    await page.getByTestId("CancelIcon").click();
    await expect(checkCellTypes(page)).rejects.toThrow();
  });
});
/**
 * *******************************************
 * Helper Functions
 * *******************************************
 */

/**
 * loadPageAndTissues
 * Load the WMG page and wait for the tissue nodes to be visible
 */
async function loadPageAndTissues(page: Page) {
  await goToWMG(page);
  await expect(page.getByTestId(TISSUE_NODE_TEST_ID)).not.toHaveCount(0);
}

/**
 * clickIntoFilter
 * Click into the filter and wait for the tooltip to be visible
 */
async function clickIntoFilter(page: Page, filterName: string) {
  await page.getByTestId(filterName).getByRole("button").first().click();
  await page.getByRole("tooltip").waitFor();
}

/**
 * filterCellType
 * Filter cell type 'B cell' and exit filter mode
 */
async function filterCellType(page: Page) {
  await page.getByRole("combobox").first().click();
  await page
    .getByRole("option", { name: CELL_TYPE_FILTER, exact: true })
    .click();
  await page.keyboard.press("Escape");
}

/**
 * filterTissues
 * Filter the tissuesfrom the left panel
 */
async function filterTissues(page: Page, filteredTissues: string[]) {
  await clickIntoFilter(page, TISSUE_FILTER_TEST_ID);
  for (const tissue of filteredTissues) {
    await page.getByRole("option", { name: tissue }).click();
  }
  await page.keyboard.press("Escape");
}

/**
 * checkTissues
 * Check that the filtered tissues are visible, and expanded tissues are expanded
 */
async function checkTissues(
  page: Page,
  filteredTissues: string[],
  expandedTissues: string[]
) {
  if (filteredTissues.length !== 0) {
    await expect(page.getByTestId("tissue-name")).toHaveCount(
      filteredTissues.length
    );
  }

  let i = 0;
  for (const tissue of expandedTissues) {
    await expect(page.getByTestId(`cell-type-labels-${tissue}`)).toBeVisible();
    const height = await page
      .getByTestId(`cell-type-labels-${tissue}`)
      .getAttribute("height");
    toInteger(height) > 20 && i++;
  }
  await expect(i).toEqual(expandedTissues.length);
}

/**
 * checkCellTypes
 * Check that all cells under expanded tissues are 'B Cell'
 */
async function checkCellTypes(page: Page) {
  const cells = await page.getByTestId("cell-type-name").allInnerTexts();
  for (const cell of cells) {
    await expect(cell).toEqual(CELL_TYPE_FILTER);
  }
}
