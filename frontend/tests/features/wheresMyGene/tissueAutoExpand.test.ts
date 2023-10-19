import { expect, Page } from "@playwright/test";
import { toInteger } from "lodash";
import { test } from "tests/common/test";
import { collapseTissue, expandTissue, tryUntil } from "tests/utils/helpers";
import {
  goToWMG,
  WMG_WITH_SEEDED_GENES_AND_CELL_TYPES,
  WMG_WITH_SEEDED_GENES_AND_TISSUES,
} from "tests/utils/wmgUtils";

const FILTERED_TISSUES = ["abdomen", "axilla", "blood"];
const TISSUE_NODE_TEST_ID = "tissue-name";
const TISSUE_FILTER_LABEL = "Tissue";
const TISSUE_FILTER_TEST_ID = "tissue-filter";
const CELL_TYPE_FILTER_TEST_ID = "celltype-filter";
const CELL_TYPE_FILTERS = ["B cell", "B-1a B cell", "B-1b B cell"];
const CELL_TYPE_TEST_ID = "cell-type-name";
const SEX_FILTER_TEST_ID = "sex-filter";
const SEX_FILTER_LABEL = "Sex";
const SEX_FILTER_SELECTION = "female";
const PUBLICATION_FILTER_TEST_ID = "publication-filter";
const PUBLICATION_FILTER_LABEL = "Publication";
const PUBLICATION_FILTER_SELECTION = [
  "Ahern et al. (2022) Cell",
  "Arutyunyan et al. (2023) Nature",
];
const SELF_REPORTED_ETHNICITY_FILTER_TEST_ID = "self-reported-ethnicity-filter";
const SELF_REPORTED_ETHNICITY_FILTER_LABEL = "Self-Reported Ethnicity";
const SELF_REPORTED_ETHNICITY_FILTER_SELECTION = "African";
const SELF_REPORTED_ETHNICITY_TISSUE = ["breast", "nose"];
const DISEASE_FILTER_TEST_ID = "disease-filter";
const DISEASE_FILTER_LABEL = "Disease";
const DISEASE_FILTER_SELECTION = "influenza";
const DATASET_FILTER_TEST_ID = "dataset-filter";
const DATASET_FILTER_LABEL = "Dataset";
const DATASET_FILTER_SELECTION =
  "Combined samples HTAN MSK - Single cell profiling reveals novel tumor and myeloid subpopulations in small cell lung cancer";
const DATASET_FILTER_SELECTED = "Combined samples";
const DATASET_FILTER_FILTERED_TISSUES = ["axilla", "brain"];
const EXPECTED_FILTERED_TISSUES_WITH_SEX_FILTER = ["abdomen", "blood"];
const EXPECTED_EXPANDED_TISSUES = ["blood"];
const EXPECTED_VISIBLE_CELL = ["B Cell"];
const EXPECTED_FILTERED_TISSUES_WITH_DATASET_FILTER = ["axilla", "brain"];
const EXPECTED_FILTERED_TISSUES_WITH_DISEASE_FILTER = ["blood"];
const EXPECTED_FILTERED_TISSUES_WITH_SHARE_LINK = ["lung"];

const WAIT_FOR_REQUEST_TIMEOUT_MS = 30 * 1000; // 30 seconds

const { describe } = test;

describe("WMG tissue auto-expand", () => {
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
    await filterCellTypes(page, CELL_TYPE_FILTERS.slice(0, 1));
    await checkTissues(page);
    await checkElementVisible(page, EXPECTED_VISIBLE_CELL, CELL_TYPE_TEST_ID);
    await removeCellFilter(page);
    await expect(
      checkElementVisible(page, CELL_TYPE_FILTERS, CELL_TYPE_TEST_ID)
    ).rejects.toThrow();
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
    await loadPageAndTissues(page);
    await filterTissues(page, tissues);
    await collapseTissue(page, tissues[0]);
    await collapseTissue(page, tissues[1]);
    await filterCellTypes(page, CELL_TYPE_FILTERS.slice(0, 1));
    await checkTissues(page, tissues);
    await checkElementVisible(page, EXPECTED_VISIBLE_CELL, CELL_TYPE_TEST_ID);
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
    await filterCellTypes(page, CELL_TYPE_FILTERS);
    await checkTissues(page);
    await removeCellFilter(page);
    await checkTissues(page, ["lung"], []);
    await checkElementVisible(page, cells, CELL_TYPE_TEST_ID);
  });

  /**
   * Tissue auto expansion - cross filter with Sex filter selected
   * Filter Female from Sex filter. Filter Abdomen and Blood from the Tissue
   * filter - those should be expanded, Add Cell type filtering,
   * only B Cell type should appear and expanded.
   */
  test("Tissue auto expansion - cross filter with Sex filter selected", async ({
    page,
  }) => {
    await loadPageAndTissues(page);
    await filterSelection(
      page,
      SEX_FILTER_TEST_ID,
      SEX_FILTER_LABEL,
      SEX_FILTER_SELECTION
    );
    await filterTissues(page, EXPECTED_FILTERED_TISSUES_WITH_SEX_FILTER);
    await filterCellTypes(page, CELL_TYPE_FILTERS.slice(0, 1));
    await checkTissues(page, EXPECTED_FILTERED_TISSUES_WITH_SEX_FILTER);
    await checkElementVisible(page, EXPECTED_VISIBLE_CELL, CELL_TYPE_TEST_ID);
  });

  /**
   * Tissue auto expansion - cross filter with Sex filter, check expansion
   * Filter 3 Tissues ["abdomen", "axilla", "blood"]. Collapse Abdomen. Select
   * Female from the Sex filter. Tissue filter should now only have Abdomen and
   * Blood selected. Only Abdomen and Blood should be visible. Abdomen should
   * remain collapsed. Remove Sex filter. Tissue filter should now only have
   * Abdomen and Blood selected. Only Abdomen and Blood should be visible.
   * Abdomen should remain collapsed.
   */
  test("Tissue auto expansion - cross filter with Sex filter, check expansion", async ({
    page,
  }) => {
    await loadPageAndTissues(page);
    await filterTissues(page, FILTERED_TISSUES);
    await collapseTissue(page, FILTERED_TISSUES[0]);
    await filterSelection(
      page,
      SEX_FILTER_TEST_ID,
      SEX_FILTER_LABEL,
      SEX_FILTER_SELECTION
    );
    await expect(page.getByTestId(TISSUE_NODE_TEST_ID)).toHaveCount(
      EXPECTED_FILTERED_TISSUES_WITH_SEX_FILTER.length
    );
    await checkTissues(
      page,
      EXPECTED_FILTERED_TISSUES_WITH_SEX_FILTER,
      EXPECTED_EXPANDED_TISSUES
    );
  });

  /**
   * Tissue auto expansion - cross filter with Publication filter selected
   * Filter 'Ahren et al. Cell 2022' from the Publication filter. Filter Blood
   * from the tissue filter. Blood should be expanded. Collapse Blood. add Cell
   * filter for B Cell. Only B Cell should appear and Blood should be expanded.
   */
  test("Tissue auto expansion - cross filter with Publication filter selected", async ({
    page,
  }) => {
    const tissues = ["blood"];
    await loadPageAndTissues(page);
    await filterSelection(
      page,
      PUBLICATION_FILTER_TEST_ID,
      PUBLICATION_FILTER_LABEL,
      PUBLICATION_FILTER_SELECTION[0]
    );
    await filterTissues(page, tissues);
    await collapseTissue(page, tissues[0]);
    await filterCellTypes(page, CELL_TYPE_FILTERS.slice(0, 1));
    await checkTissues(page, tissues);
    await checkElementVisible(page, EXPECTED_VISIBLE_CELL, CELL_TYPE_TEST_ID);
  });

  /**
   * Tissue auto expansion - cross filter with Publication filter, check expansion
   * Filter 3 Tissues ["abdomen", "axilla", "blood"]. Collapse Blood. Filter
   * 'Ahren et al. Neuron 2021' from the Publication filter. Only Blood should remain
   * visible and collapsed. Add Cell filter for B Cell. Only B Cell should be
   * visible under expanded Blood tissue node
   */
  test("Tissue auto expansion - cross filter with Publication filter, check expansion", async ({
    page,
  }) => {
    await loadPageAndTissues(page);
    await filterTissues(page);
    await collapseTissue(page, "blood");
    await filterSelection(
      page,
      PUBLICATION_FILTER_TEST_ID,
      PUBLICATION_FILTER_LABEL,
      PUBLICATION_FILTER_SELECTION[1]
    );
    await expect(page.getByTestId(TISSUE_NODE_TEST_ID)).toHaveCount(1);
    await checkTissues(page, ["blood"], []);
  });

  /**
   * Tissue auto expansion - cross filter with Self-Reported Ethnicity filter
   * Expand Breast. Select African from Ethnicity filter. Breast should remain expanded.
   * Add Breast and Nose to the Tissue filter. Both should appear expanded.
   */
  test("Tissue auto expansion - cross filter with Self-Reported Ethnicity filter", async ({
    page,
  }) => {
    await loadPageAndTissues(page);
    await expandTissue(page, SELF_REPORTED_ETHNICITY_TISSUE[0]);
    await filterSelection(
      page,
      SELF_REPORTED_ETHNICITY_FILTER_TEST_ID,
      SELF_REPORTED_ETHNICITY_FILTER_LABEL,
      SELF_REPORTED_ETHNICITY_FILTER_SELECTION
    );
    await filterTissues(page, SELF_REPORTED_ETHNICITY_TISSUE);
    await checkTissues(page, SELF_REPORTED_ETHNICITY_TISSUE);
  });

  /**
   * Tissue auto expansion - cross filter with Disease filter
   * Filter 3 Tissues ["abdomen", "axilla", "brain"]. From the Disease filter select influenza.
   * Blood should appear only and expanded. Add cell type filter for B Cell. Only B Cell should
   * appear under Blood and expanded. Remove influenza. Blood should remain expanded and cell count
   * increase.
   */
  test("Tissue auto expansion - cross filter with Disease filter", async ({
    page,
  }) => {
    await loadPageAndTissues(page);
    await filterTissues(page);
    await filterSelection(
      page,
      DISEASE_FILTER_TEST_ID,
      DISEASE_FILTER_LABEL,
      DISEASE_FILTER_SELECTION
    );
    await checkTissues(page, EXPECTED_FILTERED_TISSUES_WITH_DISEASE_FILTER);
    await filterCellTypes(page, CELL_TYPE_FILTERS.slice(0, 1));
    await checkElementVisible(page, EXPECTED_VISIBLE_CELL, CELL_TYPE_TEST_ID);
    await checkTissues(page, EXPECTED_FILTERED_TISSUES_WITH_DISEASE_FILTER);
    await removeCellFilter(page);
    await checkTissues(page, EXPECTED_FILTERED_TISSUES_WITH_DISEASE_FILTER);
  });

  /**
   * Tissue auto expansion - cross filter with Dataset filter
   * Filter 3 Tissues ["abdomen", "axilla", "Brain"]. From the Dataset filter select Combined Samples.
   * Axilla and Brain should appear only and expanded. Add cell type filter for B Cell. Only B Cell
   * should appear under Axilla and Brain and expanded. Remove B Cell. Remove Combined Samples. Brain
   */
  test("Tissue auto expansion - cross filter with Dataset filter", async ({
    page,
  }) => {
    await loadPageAndTissues(page);
    await filterTissues(page, DATASET_FILTER_FILTERED_TISSUES);
    await filterSelection(
      page,
      DATASET_FILTER_TEST_ID,
      DATASET_FILTER_LABEL,
      DATASET_FILTER_SELECTION,
      DATASET_FILTER_SELECTED
    );
    await checkTissues(page, EXPECTED_FILTERED_TISSUES_WITH_DATASET_FILTER);
    await filterCellTypes(page, CELL_TYPE_FILTERS.slice(0, 1));
    await checkElementVisible(page, EXPECTED_VISIBLE_CELL, CELL_TYPE_TEST_ID);
    await checkTissues(page, EXPECTED_FILTERED_TISSUES_WITH_DATASET_FILTER);
    await removeCellFilter(page);
    await checkTissues(page, EXPECTED_FILTERED_TISSUES_WITH_DATASET_FILTER);
  });

  test("Share link with genes and cellTypes", async ({ page }) => {
    await Promise.all([
      /**
       * (thuang): This test asserts that the app does use the cellTypes passed
       * in the share link in a `/filters` request.
       * If this `waitForRequest` times out, it's likely because the app is NOT
       * sending a request with the cellTypes passed in the share link.
       */
      page.waitForRequest(
        (request) => {
          if (!request.url().includes("wmg/v2/filters")) return false;

          const requestBody = JSON.parse(request.postData() || "{}");

          const requestCellTypeIds = JSON.stringify(
            requestBody.filter.cell_type_ontology_term_ids
          );

          return (
            requestCellTypeIds ===
            JSON.stringify(WMG_WITH_SEEDED_GENES_AND_CELL_TYPES.cellTypeIds)
          );
        },
        { timeout: WAIT_FOR_REQUEST_TIMEOUT_MS }
      ),
      loadPageAndTissues(page, WMG_WITH_SEEDED_GENES_AND_CELL_TYPES.URL),
    ]);

    await checkElementVisible(
      page,
      WMG_WITH_SEEDED_GENES_AND_CELL_TYPES.cellTypes,
      CELL_TYPE_TEST_ID
    );
    await checkTissues(page, EXPECTED_FILTERED_TISSUES_WITH_SHARE_LINK);
  });

  test("Share link with genes and tissues", async ({ page }) => {
    await Promise.all([
      /**
       * (thuang): This test asserts that the app does use the tissues passed
       * in the share link in a `/filters` request.
       * If this `waitForRequest` times out, it's likely because the app is NOT
       * sending a request with the tissues passed in the share link.
       */
      page.waitForRequest(
        (request) => {
          if (!request.url().includes("wmg/v2/filters")) return false;

          const requestBody = JSON.parse(request.postData() || "{}");

          const requestTissueIds = JSON.stringify(
            requestBody.filter.tissue_ontology_term_ids
          );

          return (
            requestTissueIds ===
            JSON.stringify(WMG_WITH_SEEDED_GENES_AND_TISSUES.tissueIds)
          );
        },
        { timeout: WAIT_FOR_REQUEST_TIMEOUT_MS }
      ),
      loadPageAndTissues(page, WMG_WITH_SEEDED_GENES_AND_TISSUES.URL),
    ]);

    await checkTissues(page, WMG_WITH_SEEDED_GENES_AND_TISSUES.tissues);
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
async function loadPageAndTissues(page: Page, url?: string) {
  await goToWMG(page, url);
  await expect(page.getByTestId(TISSUE_NODE_TEST_ID)).not.toHaveCount(0);
}

/**
 * clickIntoFilter
 * Click into the filter and wait for the tooltip to be visible
 */
async function clickIntoFilter(
  page: Page,
  filterName: string,
  filterLabel = TISSUE_FILTER_LABEL
) {
  await page
    .getByTestId(filterName)
    .getByRole("button", { name: filterLabel, exact: true })
    .click();
  await page.getByRole("tooltip").waitFor();
}

/**
 * filterCellTypes
 * Add cell types to the cell type filter
 */
async function filterCellTypes(page: Page, cellTypes: string[]) {
  await page
    .getByTestId(CELL_TYPE_FILTER_TEST_ID)
    .getByRole("combobox")
    .click();

  for (const cellType of cellTypes) {
    await page.getByRole("option", { name: cellType, exact: true }).click();
  }

  await page.keyboard.press("Escape");
}

/**
 * filterTissues
 * Filter the tissues from the left panel
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
  await tryUntil(
    async () => {
      if (filteredTissues.length !== 0) {
        await expect(page.getByTestId(TISSUE_NODE_TEST_ID)).toHaveCount(
          filteredTissues.length
        );
        await checkElementVisible(page, filteredTissues, TISSUE_NODE_TEST_ID);
      }

      let countExpanded = 0;
      for (const tissue of expandedTissues) {
        await expect(
          page.getByTestId(`cell-type-labels-${tissue}`)
        ).toBeVisible();
        const height = await page
          .getByTestId(`cell-type-labels-${tissue}`)
          .getAttribute("height");
        toInteger(height) > 20 && countExpanded++;
      }
      await expect(countExpanded).toEqual(expandedTissues.length);
    },
    { page }
  );
}

/**
 * checkElementVisible
 * Check that all elements are visible
 */
async function checkElementVisible(
  page: Page,
  filteredElements: string[],
  testId: string
) {
  await tryUntil(
    async () => {
      const elements = await (
        await page.getByTestId(testId).allInnerTexts()
      ).map((str) => str.toLowerCase());
      expect(elements).toEqual(
        expect.arrayContaining(filteredElements.map((str) => str.toLowerCase()))
      );
    },
    { page }
  );
}

/**
 * removeCellFilter
 * Specifically remove the first cell filter (B cell)
 */
async function removeCellFilter(page: Page) {
  await page
    .getByTestId(`cell-type-tag-${CELL_TYPE_FILTERS[0]}`)
    .getByTestId("ClearIcon")
    .click();
}

/**
 * filterSelection
 * Filter the selection from the dropdown
 */
async function filterSelection(
  page: Page,
  filterTestId: string,
  filterLabel: string,
  selection: string,
  filterSelected: string = selection
) {
  await tryUntil(
    async () => {
      await clickIntoFilter(page, filterTestId, filterLabel);
      const dropDownOption = await page.getByRole("option", {
        name: selection,
        exact: true,
      });
      expect(dropDownOption).toBeVisible();
      await dropDownOption.click();
      await page.keyboard.press("Escape");
      await expect(
        page
          .getByTestId(filterTestId)
          .getByRole("button", { name: filterSelected })
      ).toBeVisible();
    },
    { page }
  );
}
