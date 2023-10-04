import { expect, test, Page } from "@playwright/test";
import { indexOf, toInteger } from "lodash";
import { collapseTissue, expandTissue } from "tests/utils/helpers";
import { conditionallyRunTests, goToWMG } from "tests/utils/wmgUtils";

const FILTERED_TISSUES = ["abdomen", "axilla", "blood"];
const TISSUE_NODE_TEST_ID = "tissue-name";
const TISSUE_FILTER_LABEL = "Tissue";
const TISSUE_FILTER_TEST_ID = "tissue-filter";
const CELL_TYPE_FILTER_TEST_ID = "celltype-filter";
const CELL_TYPE_FILTERS = ["B cell", "B-1a B cell", "B-1b B cell"];
const CELL_TYPE_TEST_ID = "cell-type-name";

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
    await filterTissues(page);
    await checkTissues(page);
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
    await filterTissues(page);
    await checkTissues(page);
    await filterTissues(page);
    await checkTissues(page, [], []);
  });
  /**
   * Filter two tissues, expect tissues to be expanded. Expect tissue list
   * to only show tissues that contain "B cell" and automatically expanded.
   * Check that only 'B Cell' cells are visible under expanded tissues.
   * When cell type filter is removed, expect all cells to be shown under
   * expanded tissues.
   */
  test("Filter cell type auto expansion - filter cell type 'B cell'", async ({
    page,
  }) => {
    await loadPageAndTissues(page);
    await filterTissues(page);
    await filterCellType(page, 1);
    await checkTissues(page);
    await checkCellTypes(page);
    await removeCellFilter(page);
    await expect(checkCellTypes(page)).rejects.toThrow();
  });
  /**
   * Filter cell type auto expansion - override tissue filter collapse state
   * Filter first 2 tissues in left panel, collapse both tissues, filter "B
   * cell". Expect both tissues to expand and only has "B cell" row
   */
  test("Filter cell type auto expansion - override tissue filter collapse state", async ({
    page,
  }) => {
    const tissues = FILTERED_TISSUES.slice(0, 2);
    const cells = CELL_TYPE_FILTERS.slice(0, 1);
    await loadPageAndTissues(page);
    await filterTissues(page, tissues);
    await collapseTissue(page, tissues[0]);
    await collapseTissue(page, tissues[1]);
    await filterCellType(page, 1);
    await checkTissues(page, tissues);
    await checkCellTypes(page, cells);
  });
  /**
   * Filter 3 tissues, filter top 3 cell types
   * (B cell, B-1a B cell, and B-1b B cell). Expect both tissues expanded.
   * Remove “B cell” from filter <-- this removes tissue filter, since only
   * “Lung” tissue has B-1a B cell and B-1b B cell
   */
  test("Filter cell type auto expansion - remove cell type filter", async ({
    page,
  }) => {
    const cells = CELL_TYPE_FILTERS.slice(1, 3);
    await loadPageAndTissues(page);
    await filterTissues(page);
    await filterCellType(page);
    await checkTissues(page);
    await removeCellFilter(page);
    await checkTissues(page, ["lung"], []);
    await checkCellTypes(page, cells);
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
  await page
    .getByTestId(filterName)
    .getByRole("button", { name: TISSUE_FILTER_LABEL, exact: true })
    .click();
  await page.getByRole("tooltip").waitFor();
}

/**
 * filterCellType
 * Filter cell type 'B cell' and exit filter mode
 */
async function filterCellType(page: Page, count = 3) {
  await page
    .getByTestId(CELL_TYPE_FILTER_TEST_ID)
    .getByRole("combobox")
    .click();
  for (let i = 0; i < count; i++) {
    await page
      .getByRole("option", { name: CELL_TYPE_FILTERS[i], exact: true })
      .click();
  }
  await page.keyboard.press("Escape");
}

/**
 * filterTissues
 * Filter the tissuesfrom the left panel
 */
async function filterTissues(
  page: Page,
  filteredTissues: string[] = FILTERED_TISSUES
) {
  await clickIntoFilter(page, TISSUE_FILTER_TEST_ID);
  for (const tissue of filteredTissues) {
    await page.getByRole("option", { name: tissue, exact: true }).click();
  }
  await page.keyboard.press("Escape");
}

/**
 * checkTissues
 * Check that only filtered tissues are visible, and expanded tissues are expanded
 */
async function checkTissues(
  page: Page,
  filteredTissues: string[] = FILTERED_TISSUES,
  expandedTissues = filteredTissues
) {
  if (filteredTissues.length !== 0) {
    await expect(page.getByTestId(TISSUE_NODE_TEST_ID)).toHaveCount(
      filteredTissues.length
    );
  }

  let countExpanded = 0;
  for (const tissue of expandedTissues) {
    await expect(page.getByTestId(`cell-type-labels-${tissue}`)).toBeVisible();
    const height = await page
      .getByTestId(`cell-type-labels-${tissue}`)
      .getAttribute("height");
    toInteger(height) > 20 && countExpanded++;
  }
  await expect(countExpanded).toEqual(expandedTissues.length);
}

/**
 * checkCellTypes
 * Check that all cells under expanded tissues are 'B Cell'
 */
async function checkCellTypes(page: Page, cellTypes = CELL_TYPE_FILTERS) {
  const cells = (
    await page.getByTestId(CELL_TYPE_TEST_ID).allInnerTexts()
  ).sort((a, b) => a.localeCompare(b));
  for (const cell of cells) {
    await expect(cell).toEqual(cellTypes[indexOf(cells, cell)]);
  }
}

/**
 * removeCellFilter
 * Remove cell filter
 */
async function removeCellFilter(page: Page) {
  await await page
    .getByTestId(`cell-type-tag-${CELL_TYPE_FILTERS[0]}`)
    .getByTestId("CancelIcon")
    .click();
}
