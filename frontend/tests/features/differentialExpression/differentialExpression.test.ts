import fs from "fs";
import { Page, expect, Locator } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { goToPage, isElementVisible, tryUntil } from "tests/utils/helpers";
import { TEST_URL } from "../../common/constants";
import { test } from "tests/common/test";
import {
  DIFFERENTIAL_EXPRESSION_METHOD_INFO_TEXT,
  DIFFERENTIAL_EXPRESSION_ORGANISM_DROPDOWN,
  DIFFERENTIAL_EXPRESSION_METHOD_DROPDOWN,
  DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX,
  DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER,
  DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER,
  DIFFERENTIAL_EXPRESSION_FILTER_TAG_PRIMARY,
  DIFFERENTIAL_EXPRESSION_FILTER_TAG_GRAY,
  DIFFERENTIAL_EXPRESSION_INSTRUCTIONS_SIDEBAR,
  DIFFERENTIAL_EXPRESSION_FIND_GENES_BUTTON,
  DIFFERENTIAL_EXPRESSION_CLEAR_ALL_BUTTON,
  DIFFERENTIAL_EXPRESSION_COPY_FILTERS_BUTTON_PREFIX,
  DIFFERENTIAL_EXPRESSION_RESULTS_DOWNLOAD_BUTTON,
  DIFFERENTIAL_EXPRESSION_SOURCE_DATA_BUTTON,
  DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_INFO,
  DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_INFO,
  DIFFERENTIAL_EXPRESSION_OPEN_IN_GE_1_BUTTON,
  DIFFERENTIAL_EXPRESSION_OPEN_IN_GE_2_BUTTON,
  DIFFERENTIAL_EXPRESSION_RESULTS_TABLE,
  DIFFERENTIAL_EXPRESSION_GENES_FILTER,
  DIFFERENTIAL_EXPRESSION_LFC_FILTER,
  DIFFERENTIAL_EXPRESSION_EFFECT_SIZE_FILTER,
  DIFFERENTIAL_EXPRESSION_RESULTS_CALLOUT,
  DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT,
  DIFFERENTIAL_EXPRESSION_FILTERS_LOADING_SPINNER,
  DIFFERENTIAL_EXPRESSION_SORT_DIRECTION,
  DIFFERENTIAL_EXPRESSION_SOURCE_DATA_SIDEBAR,
} from "src/views/DifferentialExpression/common/constants";

const { describe } = test;

describe("Differential Expression", () => {
  test.beforeEach(async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.DE}`, page);
    await waitForFiltersEndpoint(page);
  });
  describe("Query Builder", () => {
    test("All components visible", async ({ page }) => {
      await isElementVisible(
        page,
        DIFFERENTIAL_EXPRESSION_INSTRUCTIONS_SIDEBAR
      );
      await isElementVisible(page, DIFFERENTIAL_EXPRESSION_METHOD_INFO_TEXT);
      await isElementVisible(page, DIFFERENTIAL_EXPRESSION_ORGANISM_DROPDOWN);
      await isElementVisible(page, DIFFERENTIAL_EXPRESSION_METHOD_DROPDOWN);
      await isElementVisible(page, DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER);
      await isElementVisible(page, DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER);
    });

    test("Using a filter does not crossfilter into itself", async ({
      page,
    }) => {
      const filterAutocomplete = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );
      await openAutocompleteDropdown(filterAutocomplete);

      const initialDropdownItems =
        await getAutocompleteDropdownItemsCount(filterAutocomplete);

      await clickOnAutocompleteDropdownItem(filterAutocomplete, "lung");

      const finalDropdownItems =
        await getAutocompleteDropdownItemsCount(filterAutocomplete);

      expect(initialDropdownItems).toBe(finalDropdownItems);
    });
    test("Using a filter crossfilters other filters", async ({ page }) => {
      const cellTypeFilterAutocomplete = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Cell Type`
        );
      await openAutocompleteDropdown(cellTypeFilterAutocomplete);

      const initialCellTypeDropdownItems =
        await getAutocompleteDropdownItemsCount(cellTypeFilterAutocomplete);

      const tissueFilterAutocomplete = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );
      await clickOnAutocompleteDropdownItem(tissueFilterAutocomplete, "lung");

      await openAutocompleteDropdown(cellTypeFilterAutocomplete);

      const finalCellTypeDropdownItems =
        await getAutocompleteDropdownItemsCount(cellTypeFilterAutocomplete);

      expect(finalCellTypeDropdownItems).toBeLessThan(
        initialCellTypeDropdownItems
      );
    });

    test("Using a filter updates the number of cells", async ({ page }) => {
      // Get initial cell counts
      const initialCellCountGroup1 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();
      const initialCellCountGroup2 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();

      // Input "lung" for the tissue filter on group 1
      const tissueFilterAutocompleteGroup1 = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );
      await clickOnAutocompleteDropdownItem(
        tissueFilterAutocompleteGroup1,
        "lung"
      );

      // Get cell counts after applying the tissue filter
      const finalCellCountGroup1 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();
      const finalCellCountGroup2 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();

      // Assert that group 1 cell count changed and group 2 cell count remained unchanged
      expect(initialCellCountGroup1).not.toBe(finalCellCountGroup1);
      expect(initialCellCountGroup2).toBe(finalCellCountGroup2);

      // Input "neuron" in the cell type filter on group 2
      const cellGroup2FilterAutocomplete = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Cell Type`
        );
      await clickOnAutocompleteDropdownItem(
        cellGroup2FilterAutocomplete,
        "neuron"
      );

      // Get the final cell count for group 2
      const neuronCellCountGroup2 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();

      // Assert that group 2 cell count changed after applying the neuron filter
      expect(finalCellCountGroup2).not.toBe(neuronCellCountGroup2);
    });

    test("Clearing a tag from filter unselects the filter", async ({
      page,
    }) => {
      // Record initial cell count in group 1
      const initialCellCountGroup1 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();

      // Select "lung" in tissue filter for cell group 1
      const tissueFilterAutocompleteGroup1 = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );
      await clickOnAutocompleteDropdownItem(
        tissueFilterAutocompleteGroup1,
        "lung"
      );

      // Record cell count in group 1 after applying the filter
      const lungCellCountGroup1 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();

      // Click the primary tag with "lung" to remove the filter
      await tissueFilterAutocompleteGroup1
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_TAG_PRIMARY)
        .locator("span")
        .click(); // no need to wait for filters endpoint because it's cached

      // Record cell count in group 1 after removing the filter
      const finalCellCountGroup1 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();

      // Assert that cell counts change and revert after adding/deleting the filter
      expect(initialCellCountGroup1).not.toBe(lungCellCountGroup1);
      expect(lungCellCountGroup1).not.toBe(finalCellCountGroup1);
      expect(initialCellCountGroup1).toBe(finalCellCountGroup1);
    });

    test("Clearing the entire filter with the clear indicator works", async ({
      page,
    }) => {
      // Record initial cell count in group 1
      const initialCellCountGroup1 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();

      // Select "lung" in tissue filter for cell group 1
      const tissueFilterAutocompleteGroup1 = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );
      await clickOnAutocompleteDropdownItem(
        tissueFilterAutocompleteGroup1,
        "lung"
      );

      // Record cell count in group 1 after applying the filter
      const lungCellCountGroup1 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();

      // Click the ClearIndicator to remove the filter
      // no need to wait for filters endpoint because it's cached
      await tissueFilterAutocompleteGroup1.getByTestId("CloseIcon").click();

      // Record cell count in group 1 after removing the filter
      const finalCellCountGroup1 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();

      // Assert that cell counts change and revert after adding/deleting the filter
      expect(initialCellCountGroup1).not.toBe(lungCellCountGroup1);
      expect(lungCellCountGroup1).not.toBe(finalCellCountGroup1);
      expect(initialCellCountGroup1).toBe(finalCellCountGroup1);
    });

    test("Copy button copies filter over to Cell Group 2, and cell group 2 number of cells updates", async ({
      page,
    }) => {
      // Record initial cell count in group 2
      const initialCellCountGroup2 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();

      // Select "lung" in tissue filter for cell group 1
      const tissueFilterAutocompleteGroup1 = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );
      await clickOnAutocompleteDropdownItem(
        tissueFilterAutocompleteGroup1,
        "lung"
      );

      // Click the Copy button to copy the filter to Cell Group 2
      const copyButton = page.getByTestId(
        `${DIFFERENTIAL_EXPRESSION_COPY_FILTERS_BUTTON_PREFIX}tissues`
      );
      await copyButton.click();

      // Ensure the "lung" tag is present in Cell Group 2 tissue filter
      const tissueFilterGroup2 = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );

      expect(
        tissueFilterGroup2.getByTestId(
          DIFFERENTIAL_EXPRESSION_FILTER_TAG_PRIMARY
        )
      ).toBeVisible();

      // Record cell count in group 2 after applying the filter
      const lungCellCountGroup2 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();

      // Ensure cell count in group 2 updates
      expect(initialCellCountGroup2).not.toBe(lungCellCountGroup2);
    });

    test("Crossfiltering a tag that is no longer valid grays it out, does not remove it", async ({
      page,
    }) => {
      // Add "axilla" and "brain" to the tissue filter
      const tissueFilterAutocompleteGroup1 = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );
      await clickOnAutocompleteDropdownItem(tissueFilterAutocompleteGroup1, [
        "axilla",
        "brain",
      ]);

      // Add "astrocyte" to the cell type filter
      const cellTypeFilterAutocompleteGroup1 = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Cell Type`
        );
      await clickOnAutocompleteDropdownItem(
        cellTypeFilterAutocompleteGroup1,
        "astrocyte"
      );

      // Ensure "axilla" is now a gray tag (not primary)
      const axillaTag = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_TAG_GRAY)
        .filter({ hasText: "axilla" });

      expect(axillaTag).toBeVisible();
    });

    test("Switching organism clears filter selections and updates number of cells", async ({
      page,
    }) => {
      // Select "lung" in the tissue filter for group 1
      const tissueFilterAutocompleteGroup1 = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );
      await clickOnAutocompleteDropdownItem(
        tissueFilterAutocompleteGroup1,
        "lung"
      );

      // Get initial cell counts after selecting "lung"
      const initialCellCountGroup1 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();

      // Switch organism dropdown to "Mus musculus"
      const organismDropdown = page.getByTestId(
        DIFFERENTIAL_EXPRESSION_ORGANISM_DROPDOWN
      );
      await organismDropdown.click();
      await page.keyboard.press("ArrowDown");
      await page.keyboard.press("Enter");
      await waitForFiltersEndpoint(page);

      // Ensure "lung" tag is no longer present
      const lungTag = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_TAG_PRIMARY)
        .filter({ hasText: "lung" });
      await expect(lungTag).toHaveCount(0);

      // Ensure cell counts changed
      const updatedCellCountGroup1 = await page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
        .textContent();
      expect(initialCellCountGroup1).not.toBe(updatedCellCountGroup1);
    });

    test("Clear all button clears query groups", async ({ page }) => {
      // Populate "lung" in cell group 1 filter tissue
      const tissueFilterAutocompleteGroup1 = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );
      await clickOnAutocompleteDropdownItem(
        tissueFilterAutocompleteGroup1,
        "lung"
      );

      // Populate "blood" in cell group 2 filter tissue
      const tissueFilterAutocompleteGroup2 = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );
      await clickOnAutocompleteDropdownItem(
        tissueFilterAutocompleteGroup2,
        "blood"
      );

      // Ensure "Find Genes" button is enabled
      const findGenesButton = page.getByTestId(
        DIFFERENTIAL_EXPRESSION_FIND_GENES_BUTTON
      );
      await expect(findGenesButton).toBeEnabled();

      // Click "Clear all" button
      const clearAllButton = page.getByTestId(
        DIFFERENTIAL_EXPRESSION_CLEAR_ALL_BUTTON
      );
      await clearAllButton.click();
      await waitForFiltersEndpoint(page);

      // Ensure "Find Genes" button is disabled again
      await expect(findGenesButton).toBeDisabled();

      // Ensure "lung" tag is no longer present in cell group 1 filter
      const lungTag = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_TAG_PRIMARY)
        .filter({ hasText: "lung" });
      await expect(lungTag).toHaveCount(0);

      // Ensure "blood" tag is no longer present in cell group 2 filter
      const bloodTag = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER)
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_TAG_PRIMARY)
        .filter({ hasText: "blood" });
      await expect(bloodTag).toHaveCount(0);
    });
    test("Find genes button only active when both filters populated", async ({
      page,
    }) => {
      // Ensure "Find Genes" button is disabled initially
      const findGenesButton = page.getByTestId(
        DIFFERENTIAL_EXPRESSION_FIND_GENES_BUTTON
      );
      await expect(findGenesButton).toBeDisabled();

      // Populate "lung" in cell group 1 filter tissue
      const tissueFilterAutocompleteGroup1 = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );
      await clickOnAutocompleteDropdownItem(
        tissueFilterAutocompleteGroup1,
        "lung"
      );

      // Ensure "Find Genes" button is still disabled
      await expect(findGenesButton).toBeDisabled();

      // Populate "blood" in cell group 2 filter tissue
      const tissueFilterAutocompleteGroup2 = page
        .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER)
        .getByTestId(
          `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`
        );
      await clickOnAutocompleteDropdownItem(
        tissueFilterAutocompleteGroup2,
        "blood"
      );

      // Ensure "Find Genes" button is enabled
      await expect(findGenesButton).toBeEnabled();

      // Delete "lung" in cell group 1 filter tissue
      await tissueFilterAutocompleteGroup1
        .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_TAG_PRIMARY)
        .locator("span")
        .click();
      await waitForFiltersEndpoint(page);
      // Ensure "Find Genes" button is disabled again
      await expect(findGenesButton).toBeDisabled();
    });
  });

  describe("Results", () => {
    test("All tests", async ({ page }) => {
      await runDEQuery(page);

      await test.step("Cell Group 1 and 2 contain the correct number of cells and filter tags", async () => {
        // Check number of cells
        const cellGroup1CellCount = (
          await page
            .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_INFO)
            .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
            .textContent()
        )
          ?.split(" |")
          .at(0);
        expect(cellGroup1CellCount).toBeDefined();
        const filterGroup1CellCount = await page
          .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
          .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
          .textContent();
        expect(cellGroup1CellCount).toBe(filterGroup1CellCount);

        const cellGroup2CellCount = (
          await page
            .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_INFO)
            .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
            .textContent()
        )
          ?.split(" |")
          .at(0);
        expect(cellGroup2CellCount).toBeDefined();
        const filterGroup2CellCount = await page
          .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER)
          .getByTestId(DIFFERENTIAL_EXPRESSION_FILTER_CELL_COUNT)
          .textContent();
        expect(cellGroup2CellCount).toBe(filterGroup2CellCount);

        const cellGroup1Chips = await page
          .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_INFO)
          .locator('[class*="MuiChip-label"]')
          .all();
        expect(cellGroup1Chips).toHaveLength(2);
        expect(cellGroup1Chips[0]).toHaveText("1 tissue");
        expect(cellGroup1Chips[1]).toHaveText("1 cell type");

        const cellGroup2Chips = await page
          .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_INFO)
          .locator('[class*="MuiChip-label"]')
          .all();
        expect(cellGroup2Chips).toHaveLength(2);
        expect(cellGroup2Chips[0]).toHaveText("1 tissue");
        expect(cellGroup2Chips[1]).toHaveText("1 cell type");
      });

      await test.step("Open in GE opens in a new tab with expected URL for Cell Group 1", async () => {
        // Click on GE button for cell group 1
        const [newPage1] = await Promise.all([
          page.waitForEvent("popup"),
          page.getByTestId(DIFFERENTIAL_EXPRESSION_OPEN_IN_GE_1_BUTTON).click(),
        ]);

        const newPageUrl1 = newPage1.url();
        expect(newPageUrl1).toContain("MZB1");
        expect(newPageUrl1).toContain("tissues=UBERON%3A0002048");
        expect(newPageUrl1).toContain("cellTypes=plasma+cell");
        expect(newPageUrl1).toContain("ver=2");
        await newPage1.close();

        // Click on GE button for cell group 2
        const [newPage2] = await Promise.all([
          page.waitForEvent("popup"),
          page.getByTestId(DIFFERENTIAL_EXPRESSION_OPEN_IN_GE_2_BUTTON).click(),
        ]);

        const newPageUrl2 = newPage2.url();
        expect(newPageUrl2).toContain("SAA1");
        expect(newPageUrl2).toContain("tissues=UBERON%3A0002048");
        expect(newPageUrl2).toContain("cellTypes=acinar+cell");
        expect(newPageUrl2).toContain("ver=2");
        await newPage2.close();
      });

      // TODO: Figure out how to test when callout is present
      await test.step("No overlapping cells info callout present", async () => {
        await expect(
          page.getByTestId(DIFFERENTIAL_EXPRESSION_RESULTS_CALLOUT)
        ).not.toBeVisible();
      });

      await test.step("Table filters work properly", async () => {
        const lastPage = page
          .locator(".MuiPagination-ul")
          .locator("li")
          .nth(-2);
        const pageCountBefore = await lastPage.textContent();
        const effectSizeFilterInput = page
          .getByTestId(DIFFERENTIAL_EXPRESSION_EFFECT_SIZE_FILTER)
          .locator("input");
        await effectSizeFilterInput.fill(">1.0");

        const pageCountAfterEffectSize = await lastPage.textContent();
        expect(parseFloat(pageCountAfterEffectSize ?? "0")).toBeLessThan(
          parseFloat(pageCountBefore ?? "0")
        );
        await effectSizeFilterInput.fill("");

        const lfcFilterInput = page
          .getByTestId(DIFFERENTIAL_EXPRESSION_LFC_FILTER)
          .locator("input");
        await lfcFilterInput.fill("<-.1");

        const pageCountAfterLfc = await lastPage.textContent();
        expect(parseFloat(pageCountAfterLfc ?? "0")).toBeLessThan(
          parseFloat(pageCountBefore ?? "0")
        );
        await lfcFilterInput.fill("");

        expect(parseFloat(pageCountAfterLfc ?? "0")).toBeLessThan(
          parseFloat(pageCountBefore ?? "0")
        );

        const geneFilterInput = page
          .getByTestId(DIFFERENTIAL_EXPRESSION_GENES_FILTER)
          .locator("input");
        await geneFilterInput.fill("JCHAIN");

        const rowCountAfterGeneFilter = await page
          .getByTestId(DIFFERENTIAL_EXPRESSION_RESULTS_TABLE)
          .locator("tbody tr")
          .count();
        expect(rowCountAfterGeneFilter).toBe(1);
        await geneFilterInput.fill("");
      });

      await test.step("Sorting by effect size works", async () => {
        const effectSizeHeader = page.getByTestId(
          DIFFERENTIAL_EXPRESSION_SORT_DIRECTION
        );
        const effectSizeColumn = page
          .getByTestId(DIFFERENTIAL_EXPRESSION_RESULTS_TABLE)
          .locator("tbody tr td:nth-child(3)");
        const firstEffectSizeValue = await effectSizeColumn
          .first()
          .textContent();
        expect(parseFloat(firstEffectSizeValue ?? "0")).toBeGreaterThan(0);

        // Click to sort by effect size
        await effectSizeHeader.click();
        await tryUntil(
          async () => {
            // Get the first value in the effect size column after sorting
            const firstEffectSizeValueAfterFirstClick = await effectSizeColumn
              .first()
              .textContent();
            expect(
              parseFloat(firstEffectSizeValueAfterFirstClick ?? "0")
            ).toBeLessThan(0);
          },
          { page }
        );

        // Click again to sort by effect size in the opposite direction
        await effectSizeHeader.click();
        await tryUntil(
          async () => {
            const firstEffectSizeValueAfterSecondClick = await effectSizeColumn
              .first()
              .textContent();
            expect(
              parseFloat(firstEffectSizeValueAfterSecondClick ?? "0")
            ).toBeGreaterThan(0);
          },
          { page }
        );
        // Get the first value in the effect size column after sorting again
      });

      await test.step("Source Data button reveals sidebar with expected sections and at least one collection each", async () => {
        // Click the Source Data button
        const sourceDataButton = page.getByTestId(
          DIFFERENTIAL_EXPRESSION_SOURCE_DATA_BUTTON
        );
        await sourceDataButton.click();

        // Assert that the sidebar is visible
        const sidebar = page.getByTestId(
          DIFFERENTIAL_EXPRESSION_SOURCE_DATA_SIDEBAR
        );
        await expect(sidebar).toBeVisible();

        // Assert that "Cell Group 1 Collections" text is visible
        const cellGroup1Collections = sidebar.getByText(
          "Cell Group 1 Collections"
        );
        await expect(cellGroup1Collections).toBeVisible();
        const cellGroup1CollectionsHeader =
          cellGroup1Collections.locator("+ div");

        await expect(cellGroup1CollectionsHeader).toBeVisible();

        const cellGroup1CollectionsTable =
          cellGroup1CollectionsHeader.locator("+ div");

        await expect(cellGroup1CollectionsTable).toBeVisible();
        // Assert that "Cell Group 2 Collections" text is visible
        const cellGroup2Collections = sidebar.getByText(
          "Cell Group 2 Collections"
        );
        await expect(cellGroup2Collections).toBeVisible();
        const cellGroup2CollectionsHeader =
          cellGroup2Collections.locator("+ div");

        await expect(cellGroup2CollectionsHeader).toBeVisible();
        const cellGroup2CollectionsTable =
          cellGroup2CollectionsHeader.locator("+ div");

        await expect(cellGroup2CollectionsTable).toBeVisible();

        // Close the sidebar
        await page.getByText("Source Data").locator("+ button").click();
      });

      await test.step("Download button downloads the currently filtered CSV and has the expected content", async () => {
        // Filter for MZB1, JCHAIN in gene filter
        await page
          .getByTestId(DIFFERENTIAL_EXPRESSION_GENES_FILTER)
          .locator("input")
          .fill("JCHAIN,MZB1,IGKC");

        // Click the download button
        const downloadButton = page.getByTestId(
          DIFFERENTIAL_EXPRESSION_RESULTS_DOWNLOAD_BUTTON
        );
        const [download] = await Promise.all([
          page.waitForEvent("download"),
          downloadButton.click(),
        ]);

        // Read the downloaded CSV file
        const path = await download.path();
        const csvContent = fs.readFileSync(path, "utf-8");
        const csvRows = csvContent.split("\n");

        // Ensure the first two rows contain the proper filters
        expect(csvRows[0]).toEqual(
          "# Query Group 1 Filters: tissues: lung | cellTypes: plasma cell"
        );
        expect(csvRows[1]).toEqual(
          "# Query Group 2 Filters: tissues: lung | cellTypes: acinar cell"
        );

        // Ensure the next row contains column names
        const columnNames = csvRows[2].split(",");
        expect(columnNames).toEqual([
          "Gene",
          "Log Fold Change",
          "Effect Size",
          "Adjusted P-Value",
        ]);

        // Ensure the next three rows contain data
        const dataRow1 = csvRows[3].split(",");
        const dataRow2 = csvRows[4].split(",");
        const dataRow3 = csvRows[5].split(",");
        expect(dataRow1.length).toBe(columnNames.length);
        expect(dataRow2.length).toBe(columnNames.length);
        expect(dataRow3.length).toBe(columnNames.length);
      });
    });
  });
});

const waitForFiltersEndpoint = async (page: Page) => {
  await tryUntil(
    async () => {
      await expect(
        page.getByTestId(DIFFERENTIAL_EXPRESSION_FILTERS_LOADING_SPINNER)
      ).toHaveCount(0);
    },
    { page }
  );
};

const getAutocompleteDropdownItemsCount = async (autocomplete: Locator) => {
  await autocomplete.locator("~ * li").first().waitFor();
  return await autocomplete.locator("~ * li").count();
};

const openAutocompleteDropdown = async (autocomplete: Locator) => {
  await tryUntil(
    async () => {
      await autocomplete.click();
      expect(await autocomplete.locator("~ * li").count()).toBeGreaterThan(0);
    },
    { page: autocomplete.page() }
  );
};

const clickOnAutocompleteDropdownItem = async (
  autocomplete: Locator,
  itemText: string | string[]
) => {
  if (typeof itemText === "string") {
    itemText = [itemText];
  }
  await openAutocompleteDropdown(autocomplete);
  for (const text of itemText) {
    await autocomplete.locator("input").fill(text);
    await autocomplete.locator("~ * li", { hasText: text }).first().click();
  }

  await waitForFiltersEndpoint(autocomplete.page());
};

const runDEQuery = async (page: Page) => {
  // Type "lung" in tissue filter for group 1
  const tissueFilterAutocompleteGroup1 = page
    .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
    .getByTestId(`${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Tissue`);
  await clickOnAutocompleteDropdownItem(tissueFilterAutocompleteGroup1, "lung");

  // Hit the copy button for tissue filter in group 1
  const copyButtonGroup1 = page.getByTestId(
    `${DIFFERENTIAL_EXPRESSION_COPY_FILTERS_BUTTON_PREFIX}tissues`
  );
  await copyButtonGroup1.click();

  // Type "plasma cell" in cell type filter for group 1
  const cellTypeFilterAutocompleteGroup1 = page
    .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_FILTER)
    .getByTestId(
      `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Cell Type`
    );
  await clickOnAutocompleteDropdownItem(
    cellTypeFilterAutocompleteGroup1,
    "plasma cell"
  );
  // Type "acinar cell" in cell type filter for group 2
  const cellTypeFilterAutocompleteGroup2 = page
    .getByTestId(DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_FILTER)
    .getByTestId(
      `${DIFFERENTIAL_EXPRESSION_FILTER_AUTOCOMPLETE_PREFIX}Cell Type`
    );
  await clickOnAutocompleteDropdownItem(
    cellTypeFilterAutocompleteGroup2,
    "acinar cell"
  );

  // Hit the "Find Genes" button
  const findGenesButton = page.getByTestId(
    DIFFERENTIAL_EXPRESSION_FIND_GENES_BUTTON
  );
  await findGenesButton.click();

  // Ensure the results are all visible
  await isElementVisible(page, DIFFERENTIAL_EXPRESSION_RESULTS_DOWNLOAD_BUTTON);
  await isElementVisible(page, DIFFERENTIAL_EXPRESSION_SOURCE_DATA_BUTTON);
  await isElementVisible(page, DIFFERENTIAL_EXPRESSION_CELL_GROUP_1_INFO);
  await isElementVisible(page, DIFFERENTIAL_EXPRESSION_CELL_GROUP_2_INFO);
  await isElementVisible(page, DIFFERENTIAL_EXPRESSION_OPEN_IN_GE_1_BUTTON);
  await isElementVisible(page, DIFFERENTIAL_EXPRESSION_OPEN_IN_GE_2_BUTTON);
  await isElementVisible(page, DIFFERENTIAL_EXPRESSION_RESULTS_TABLE);
  await isElementVisible(page, DIFFERENTIAL_EXPRESSION_GENES_FILTER);
  await isElementVisible(page, DIFFERENTIAL_EXPRESSION_LFC_FILTER);
  await isElementVisible(page, DIFFERENTIAL_EXPRESSION_EFFECT_SIZE_FILTER);
};
