import { ElementHandle, expect, Page, test } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import type { RawPrimaryFilterDimensionsResponse } from "src/common/queries/wheresMyGene";
import { goToPage, isDevStagingProd, tryUntil } from "tests/utils/helpers";
import { TEST_URL } from "../common/constants";
import { getTestID, getText } from "../utils/selectors";
import { TISSUE_DENY_LIST } from "./fixtures/wheresMyGene/tissueRollup";

const HOMO_SAPIENS_TERM_ID = "NCBITaxon:9606";

const GENE_LABELS_ID = "gene-labels";
const CELL_TYPE_LABELS_ID = "cell-type-name";
const ADD_TISSUE_ID = "add-tissue";
const ADD_GENE_ID = "add-gene";
const GENE_DELETE_BUTTON = "gene-delete-button";
const SOURCE_DATA_BUTTON_ID = "source-data-button";
const SOURCE_DATA_LIST_SELECTOR = `[data-test-id="source-data-list"]`;

const MUI_CHIP_ROOT = ".MuiChip-root";
const FILTERS_PANEL_NOT_FOUND = "Filters panel not found";

const { describe, skip } = test;

describe("Where's My Gene", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test("renders the getting started UI", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    // Getting Started section
    await expect(page).toHaveSelector(getText("STEP 1"));
    await expect(page).toHaveSelector(getText("Add Tissues"));

    await expect(page).toHaveSelector(getText("STEP 2"));
    await expect(page).toHaveSelector(getText("Add Genes"));

    await expect(page).toHaveSelector(getText("STEP 3"));
    await expect(page).toHaveSelector(getText("Explore Gene Expression"));

    // Filters Panel
    // (thuang): `*` is for intermediate match
    // https://playwright.dev/docs/selectors#intermediate-matches

    async function getTissueSelectorButton() {
      return page.$(getTestID(ADD_TISSUE_ID));
    }

    async function getGeneSelectorButton() {
      return page.$(getTestID(ADD_GENE_ID));
    }

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    await selectFirstOption(page);

    await clickUntilOptionsShowUp(getGeneSelectorButton, page);
    await selectFirstOption(page);

    const filtersPanel = await page.$("*css=div >> text=Filters");

    await expect(filtersPanel).toHaveSelector(getText("Dataset"));
    await expect(filtersPanel).toHaveSelector(getText("Disease"));
    await expect(filtersPanel).toHaveSelector(
      getText("Self-Reported Ethnicity")
    );
    await expect(filtersPanel).toHaveSelector(getText("Sex"));

    // Legend
    const Legend = await page.$("*css=div >> text=Gene Expression");

    await expect(Legend).toHaveSelector(getText("Gene Expression"));
    await expect(Legend).toHaveSelector(getText("Expressed in Cells (%)"));

    // Info Panel
    const InfoPanel = await page.$("*css=div >> text=Filters");

    await expect(InfoPanel).toHaveSelector(getText("Methodology"));
    await expect(InfoPanel).toHaveSelector(
      getText("After filtering cells with low coverage ")
    );
    await expect(InfoPanel).toHaveSelector(getText("Source Data"));
  });

  test("Filters and Heatmap", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    async function getTissueSelectorButton() {
      return page.$(getTestID(ADD_TISSUE_ID));
    }

    async function getGeneSelectorButton() {
      return page.$(getTestID(ADD_GENE_ID));
    }

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    await selectFirstOption(page);

    await clickUntilOptionsShowUp(getGeneSelectorButton, page);
    await selectFirstOption(page);

    await waitForHeatmapToRender(page);

    const sexSelector = await getSexSelector();

    if (!sexSelector) throw Error("No sexSelector found");

    const selectedSexesBefore = await sexSelector.$$(MUI_CHIP_ROOT);

    await expect(selectedSexesBefore.length).toBe(0);

    await clickUntilOptionsShowUp(getSexSelectorButton, page);

    await selectFirstOption(page);

    const selectedSexesAfter = await sexSelector.$$(MUI_CHIP_ROOT);

    await expect(selectedSexesAfter.length).toBe(1);

    async function getFiltersPanel() {
      return page.$(getTestID("filters-panel"));
    }

    async function getSexSelector() {
      const filtersPanel = await getFiltersPanel();

      if (!filtersPanel) {
        throw Error(FILTERS_PANEL_NOT_FOUND);
      }

      return filtersPanel.$("*css=div >> text=Sex");
    }

    async function getSexSelectorButton() {
      const filtersPanel = await getFiltersPanel();

      if (!filtersPanel) {
        throw Error(FILTERS_PANEL_NOT_FOUND);
      }

      await filtersPanel.$("*css=div >> text=Sex");
      return filtersPanel.$("*css=button >> text=Sex");
    }
  });

  test("Source Data", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    async function getSourceDataButton() {
      return page.$(getTestID(SOURCE_DATA_BUTTON_ID));
    }
    async function getTissueSelectorButton() {
      return page.$(getTestID(ADD_TISSUE_ID));
    }

    async function getGeneSelectorButton() {
      return page.$(getTestID(ADD_GENE_ID));
    }

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    await selectFirstOption(page);

    await clickUntilOptionsShowUp(getGeneSelectorButton, page);
    await selectFirstOption(page);

    await waitForHeatmapToRender(page);
    await clickUntilSidebarShowsUp(getSourceDataButton, page);
    await expect(page).toHaveSelector(
      getText(
        "After filtering cells with low coverage (less than 500 genes expressed)"
      )
    );

    await tryUntil(
      async () => {
        const sourceDataList = await page.$(SOURCE_DATA_LIST_SELECTOR);

        if (!sourceDataList) throw Error("no source data displayed");

        const sourceDataListItems = await sourceDataList?.$$(
          ".MuiListItem-root"
        );

        expect(sourceDataListItems?.length).toBeGreaterThan(0);

        await page.mouse.click(0, 0);

        async function getFiltersPanel() {
          return page.$(getTestID("filters-panel"));
        }
        async function getDatasetSelector() {
          const filtersPanel = await getFiltersPanel();

          if (!filtersPanel) {
            throw Error(FILTERS_PANEL_NOT_FOUND);
          }

          return filtersPanel.$("*css=button >> text=Dataset");
        }
        const datasetSelector = await getDatasetSelector();
        if (!datasetSelector) throw Error("No datasetSelector found");
        const selectedDatasetsBefore = await datasetSelector.$$(MUI_CHIP_ROOT);
        await expect(selectedDatasetsBefore.length).toBe(0);
        await clickUntilOptionsShowUp(getDatasetSelector, page);
        await selectFirstOption(page);
        await clickUntilSidebarShowsUp(getSourceDataButton, page);

        const sourceDataListAfter = await page.$(SOURCE_DATA_LIST_SELECTOR);

        if (!sourceDataListAfter)
          throw Error(
            "no source data displayed after selecting dataset filter"
          );

        const sourceDataListAfterItems = await sourceDataListAfter?.$$(
          ".MuiListItem-root"
        );

        expect(sourceDataListAfterItems?.length).toBeGreaterThan(0);
      },
      { page }
    );
  });

  test("Hierarchical Clustering", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    async function getTissueSelectorButton() {
      return page.$(getTestID(ADD_TISSUE_ID));
    }

    async function getGeneSelectorButton() {
      return page.$(getTestID(ADD_GENE_ID));
    }

    const TISSUE_COUNT = 1;
    const GENE_COUNT = 3;

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    await selectFirstNOptions(TISSUE_COUNT, page);

    await clickUntilOptionsShowUp(getGeneSelectorButton, page);
    await selectFirstNOptions(GENE_COUNT, page);

    const beforeGeneNames = await getNames(
      `${getTestID(GENE_LABELS_ID)} span`,
      page
    );

    const beforeCellTypeNames = await getNames(
      `${getTestID(CELL_TYPE_LABELS_ID)}`,
      page
    );

    expect(beforeGeneNames.length).toBe(GENE_COUNT);

    // (thuang): Sometimes when API response is slow, we'll not capture all the
    // cell type names, so a sanity check that we expect at least 100 names
    expect(beforeCellTypeNames.length).toBeGreaterThan(100);

    const cellTypeSortDropdown = await page.locator(
      getTestID("cell-type-sort-dropdown")
    );
    await cellTypeSortDropdown.click();
    await selectNthOption(2, page);

    const geneSortDropdown = await page.locator(
      getTestID("gene-sort-dropdown")
    );
    await geneSortDropdown.click();
    await selectNthOption(2, page);

    const afterGeneNames = await getNames(
      `${getTestID(GENE_LABELS_ID)} span`,
      page
    );

    const afterCellTypeNames = await getNames(
      `${getTestID(CELL_TYPE_LABELS_ID)}`,
      page
    );

    await tryUntil(
      async () => {
        expect(afterGeneNames.length).toBe(beforeGeneNames.length);
        expect(afterCellTypeNames.length).toBe(beforeCellTypeNames.length);

        expect(afterGeneNames).not.toEqual(beforeGeneNames);
        expect(afterCellTypeNames).not.toEqual(beforeCellTypeNames);
      },
      { page }
    );
  });

  test("delete genes", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    async function getTissueSelectorButton() {
      return page.$(getTestID(ADD_TISSUE_ID));
    }

    async function getGeneSelectorButton() {
      return page.$(getTestID(ADD_GENE_ID));
    }

    await clickUntilOptionsShowUp(getTissueSelectorButton, page);
    await selectFirstNOptions(1, page);

    await clickUntilOptionsShowUp(getGeneSelectorButton, page);
    await selectFirstNOptions(3, page);

    await waitForHeatmapToRender(page);

    const beforeGeneNames = await getNames(
      `${getTestID(GENE_LABELS_ID)} span`,
      page
    );

    await tryUntil(
      async () => {
        await page.keyboard.press("Backspace");

        // Testing single gene delete
        await page.hover(".gene-label-container");
        await page.click(getTestID(GENE_DELETE_BUTTON));

        const afterGeneNames = await getNames(
          `${getTestID(GENE_LABELS_ID)} span`,
          page
        );

        expect(afterGeneNames.length).toBe(beforeGeneNames.length - 1);
        expect(afterGeneNames).not.toEqual(beforeGeneNames);
      },
      { page }
    );
  });

  describe("tissue rollup", () => {
    test("does NOT show tissues on the deny list", async ({ page }) => {
      const [response] = await Promise.all([
        page.waitForResponse("**/primary_filter_dimensions"),
        goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page),
      ]);

      const primaryFilterDimensions =
        (await response.json()) as RawPrimaryFilterDimensionsResponse;

      const humanTissues =
        primaryFilterDimensions.tissue_terms[HOMO_SAPIENS_TERM_ID];

      const tissueIds = new Set(
        humanTissues.map((tissue) => Object.keys(tissue)[0])
      );

      const hasDeniedTissue = TISSUE_DENY_LIST.some((deniedTissueId) =>
        tissueIds.has(deniedTissueId)
      );

      expect(hasDeniedTissue).toBe(false);
    });
  });
});

async function getNames(selector: string, page: Page): Promise<string[]> {
  const geneLabelsLocator = await page.locator(selector);
  await tryUntil(
    async () => {
      const names = await geneLabelsLocator.allTextContents();
      expect(typeof names[0]).toBe("string");
    },
    { page }
  );

  return geneLabelsLocator.allTextContents();
}

async function clickUntilOptionsShowUp(
  getTarget: () => Promise<ElementHandle<SVGElement | HTMLElement> | null>,
  page: Page
) {
  await tryUntil(
    async () => {
      const target = await getTarget();

      if (!target) throw Error("no target");

      await target.click();
      const tooltip = await page.$("[role=tooltip]");

      if (!tooltip) throw Error("no tooltip");

      const options = await tooltip.$$("[role=option]");

      if (!options?.length) throw Error("no options");
    },
    { page }
  );
}
async function clickUntilSidebarShowsUp(
  getTarget: () => Promise<ElementHandle<SVGElement | HTMLElement> | null>,
  page: Page
) {
  await tryUntil(
    async () => {
      const target = await getTarget();

      if (!target) throw Error("no target");

      await target.click();
      const drawer = await page.$("[class=bp4-drawer-header]");

      if (!drawer) throw Error("no drawer");
    },
    { page }
  );
}
// (thuang): This only works when a dropdown is open
async function selectFirstOption(page: Page) {
  await selectFirstNOptions(1, page);
}

async function selectFirstNOptions(count: number, page: Page) {
  for (let i = 0; i < count; i++) {
    await page.keyboard.press("ArrowDown");
    await page.keyboard.press("Enter");
  }

  await page.keyboard.press("Escape");
}

async function selectNthOption(number: number, page: Page) {
  // (thuang): Since the first option is now active, we need to offset by 1
  const step = number - 1;

  for (let i = 0; i < step; i++) {
    await page.keyboard.press("ArrowDown");
  }

  await page.keyboard.press("Enter");
  await page.keyboard.press("Escape");
}

async function waitForHeatmapToRender(page: Page) {
  await tryUntil(
    async () => {
      const canvases = await page.$$("canvas");
      await expect(canvases.length).not.toBe(0);
    },
    { page }
  );
}
