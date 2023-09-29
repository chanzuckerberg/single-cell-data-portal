import { expect, test, Page } from "@playwright/test";
import { toInteger } from "lodash";
import { collapseTissue, expandTissue } from "tests/utils/helpers";
import { conditionallyRunTests, goToWMG } from "tests/utils/wmgUtils";

const { describe } = test;
const TISSUE_NODE_TEST_ID = "tissue-name";
const TISSUE_FILTER_TEST_ID = "tissue-filter";
const CELL_TYPE_FILTER = "B cell";

describe("WMG tissue auto-expand", () => {
  conditionallyRunTests({ forceRun: true });
  /**
   * Click the first tissue to expand, filter two tissues from left
   * panel. Expect only two tissues present and both expanded
   * automatically.
   */
  test("Filter tissue auto expansion", async ({ page }) => {
    const expandedTissues = ["abdomen", "axilla"];
    await loadPageAndTissues(page);
    await expandTissue(page, expandedTissues[0]);
    await filterTissues(page, expandedTissues);
    await checkTissues(page, expandedTissues, expandedTissues);
  });
  /**
   * Filter first two tissues from the left panel, collapse the first
   * tissue, add third tissue from the left panel. Expect the first tissue
   * to stay collapsed, and the other two tissues expanded.
   */
  test("Filter tissue auto expansion - respect user collapsed tissue state", async ({
    page,
  }) => {
    const filteredTissues = ["abdomen", "axilla", "brain"];
    await loadPageAndTissues(page);
    await filterTissues(page, filteredTissues.slice(0, 2));
    await collapseTissue(page, filteredTissues[0]);
    await filterTissues(page, filteredTissues.slice(2, 3));
    await checkTissues(page, filteredTissues, filteredTissues.slice(1, 3));
  });
  /**
   * Filter two tissues, expect tissues to be expanded. Remove both tissues
   * from filter. Expect showing all tissues in collapsed state.
   */
  test("Filter tissue auto expansion - exit filter tissues mode", async ({
    page,
  }) => {
    const filteredTissues = ["abdomen", "axilla"];
    await loadPageAndTissues(page);
    await filterTissues(page, filteredTissues);
    await checkTissues(page, filteredTissues, filteredTissues);
    await filterTissues(page, filteredTissues);
    await checkTissues(page, [], []);
  });
  /**
   * Expect tissue list to only show tissues that contain "B cell" and
   * automatically expanded.
   */
  test("Filter tissue auto expansion - filter cell type 'B cell'", async ({
    page,
  }) => {
    const filteredTissues = ["abdomen", "axilla"];
    await loadPageAndTissues(page);
    await filterTissues(page, filteredTissues);
    await filterCellType(page);
    await checkTissues(page, filteredTissues, filteredTissues);
    await checkCellTypes(page);
  });
});

async function loadPageAndTissues(page: Page) {
  await goToWMG(page);
  await expect(page.getByTestId(TISSUE_NODE_TEST_ID)).not.toHaveCount(0);
}

async function clickIntoFilter(page: Page, filterName: string) {
  await page.getByTestId(filterName).getByRole("button").first().click();
  await page.getByRole("tooltip").waitFor();
}

async function filterCellType(page: Page) {
  await page.getByRole("combobox").first().click();
  await page
    .getByRole("option", { name: CELL_TYPE_FILTER, exact: true })
    .click();
  await page.keyboard.press("Escape");
}

async function filterTissues(page: Page, expandedTissues: string[]) {
  await clickIntoFilter(page, TISSUE_FILTER_TEST_ID);
  for (const tissue of expandedTissues) {
    await page.getByRole("option", { name: tissue }).click();
  }
  await page.keyboard.press("Escape");
}

async function checkTissues(
  page: Page,
  filteredTissues: string[],
  expandedTissues: string[]
) {
  // check that only n tissue nodes are displayed
  if (filteredTissues.length !== 0) {
    await expect(page.getByTestId("tissue-name")).toHaveCount(
      filteredTissues.length
    );
  }

  let i = 0;
  // check that the filtered tissues are visible and expanded
  for (const tissue of expandedTissues) {
    await expect(page.getByTestId(`cell-type-labels-${tissue}`)).toBeVisible();
    const height = await page
      .getByTestId(`cell-type-labels-${tissue}`)
      .getAttribute("height");
    toInteger(height) > 20 && i++;
  }
  await expect(i).toEqual(expandedTissues.length);
}

async function checkCellTypes(page: Page) {
  const cells = await page.getByTestId("cell-type-name").allInnerTexts();
  for (const cell of cells) {
    await expect(cell).toEqual(CELL_TYPE_FILTER);
  }
}
