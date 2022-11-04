import { ElementHandle, expect, Page, test } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import type { RawPrimaryFilterDimensionsResponse } from "src/common/queries/wheresMyGene";
import { goToPage, isDevStagingProd, tryUntil } from "tests/utils/helpers";
import { TEST_URL } from "../common/constants";
import { getTestID, getText } from "../utils/selectors";
import { TISSUE_DENY_LIST } from "./fixtures/wheresMyGene/tissueRollup";

const HOMO_SAPIENS_TERM_ID = "NCBITaxon:9606";

const GENE_LABELS_ID = "gene-labels";
const CELL_TYPE_LABELS_ID = "cell-type-labels";
const ADD_TISSUE_ID = "add-tissue";
const ADD_GENE_ID = "add-gene";

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

    // Beta callout
    await expect(page).toHaveSelector(getText("This feature is in beta"));

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
    await expect(filtersPanel).toHaveSelector(getText("Ethnicity"));
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

    await tryUntil(
      async () => {
        const canvases = await page.$$("canvas");
        await expect(canvases.length).not.toBe(0);
      },
      { page }
    );

    const sexSelector = await getSexSelector();

    if (!sexSelector) throw Error("No sexSelector found");

    const selectedSexesBefore = await sexSelector.$$(".MuiChip-root");

    await expect(selectedSexesBefore.length).toBe(0);

    await clickUntilOptionsShowUp(getSexSelectorButton, page);

    await selectFirstOption(page);

    const selectedSexesAfter = await sexSelector.$$(".MuiChip-root");

    await expect(selectedSexesAfter.length).toBe(1);

    async function getFiltersPanel() {
      return page.$(getTestID("filters-panel"));
    }

    async function getSexSelector() {
      const filtersPanel = await getFiltersPanel();

      if (!filtersPanel) {
        throw Error("Filters panel not found");
      }

      return filtersPanel.$("*css=div >> text=Sex");
    }

    async function getSexSelectorButton() {
      const filtersPanel = await getFiltersPanel();

      if (!filtersPanel) {
        throw Error("Filters panel not found");
      }

      await filtersPanel.$("*css=div >> text=Sex");
      return filtersPanel.$("*css=button >> text=Sex");
    }
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
      `${getTestID(GENE_LABELS_ID)} text`,
      page
    );

    const beforeCellTypeNames = await getNames(
      `${getTestID(CELL_TYPE_LABELS_ID)} text`,
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
      `${getTestID(GENE_LABELS_ID)} text`,
      page
    );

    const afterCellTypeNames = await getNames(
      `${getTestID(CELL_TYPE_LABELS_ID)} text`,
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

  test("delete genes and cell types", async ({ page }) => {
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

    await tryUntil(
      async () => {
        const canvases = await page.$$("canvas");
        await expect(canvases.length).not.toBe(0);
      },
      { page }
    );

    const beforeGeneNames = await getNames(
      `${getTestID(GENE_LABELS_ID)} text`,
      page
    );
    const beforeCellTypeNames = await getNames(
      `${getTestID(CELL_TYPE_LABELS_ID)} text`,
      page
    );

    await page.click(getText(beforeGeneNames[0]));
    await page.click(getText(beforeCellTypeNames[0]));

    await tryUntil(
      async () => {
        await page.focus(getTestID(GENE_LABELS_ID));
        await page.keyboard.press("Backspace");

        const afterGeneNames = await getNames(
          `${getTestID(GENE_LABELS_ID)} text`,
          page
        );
        const afterCellTypeNames = await getNames(
          `${getTestID(CELL_TYPE_LABELS_ID)} text`,
          page
        );

        expect(afterGeneNames.length).toBe(beforeGeneNames.length - 1);

        // (thuang): Sometimes when API response is slow, we'll not capture all the
        // cell type names, so a sanity check that we expect at least 100 names
        expect(beforeCellTypeNames.length).toBeGreaterThan(100);

        // (thuang): We need to half the cellTypeName count, because it's grabbing
        // Cell Count text elements as well.
        expect(afterCellTypeNames.length / 2).toBe(
          beforeCellTypeNames.length / 2 - 1
        );

        expect(afterGeneNames).not.toEqual(beforeGeneNames);
        expect(afterCellTypeNames).not.toEqual(beforeCellTypeNames);
      },
      { page }
    );

    const RESET_CELL_TYPES_BUTTON_ID = "reset-cell-types";

    await tryUntil(
      async () => {
        await page.click(getTestID(RESET_CELL_TYPES_BUTTON_ID));

        const resetCellTypesButton = await page.$(
          getTestID(RESET_CELL_TYPES_BUTTON_ID)
        );

        expect(resetCellTypesButton).toBe(null);
      },
      { page }
    );

    await tryUntil(
      async () => {
        const afterCellTypeNames = await getNames(
          `${getTestID(CELL_TYPE_LABELS_ID)} text`,
          page
        );

        expect(afterCellTypeNames.length).toBe(beforeCellTypeNames.length);
        expect(afterCellTypeNames).toEqual(beforeCellTypeNames);
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
  for (let i = 0; i < number; i++) {
    await page.keyboard.press("ArrowDown");
  }

  await page.keyboard.press("Enter");
  await page.keyboard.press("Escape");
}
