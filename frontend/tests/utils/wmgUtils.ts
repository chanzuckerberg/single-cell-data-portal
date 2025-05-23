import { ROUTES } from "src/common/constants/routes";
import {
  SOURCE_DATA_BUTTON_ID,
  SOURCE_DATA_LIST_ID,
  TEST_URL,
} from "../common/constants";
import { expect, Page } from "@playwright/test";
import { getText } from "tests/utils/selectors";
import {
  countLocator,
  expandTissue,
  getCellTypeNames,
  selectFirstOption,
  tryUntil,
  waitForLoadingSpinnerToResolve,
} from "./helpers";
import { ADD_GENE_BTN, ADD_TISSUE_ID } from "../common/constants";
import {
  ADD_GENE_SEARCH_PLACEHOLDER_TEXT,
  CELL_TYPE_SEARCH_PLACEHOLDER_TEXT,
} from "tests/utils/geneUtils";
import { TISSUE_NAME_LABEL_CLASS_NAME } from "src/views/WheresMyGeneV2/components/HeatMap/components/YAxisChart/constants";

const WMG_SEED_GENES = ["DPM1", "TNMD", "TSPAN6"];

/**
 * (thuang): Seed app state with some genes
 */
export const WMG_WITH_SEEDED_GENES = {
  URL:
    `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}?` +
    `genes=${encodeURIComponent(WMG_SEED_GENES.join(","))}&ver=2`,
  genes: WMG_SEED_GENES,
};

const WMG_SEED_CELL_TYPES = ["B-1a B cell"];
/**
 * B-1a B cell's ontology id is: CL:0000820
 */
const WMG_SEED_CELL_TYPE_IDS = ["CL:0000820"];

export const WMG_WITH_SEEDED_GENES_AND_CELL_TYPES = {
  URL:
    `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}?` +
    `genes=${encodeURIComponent(WMG_SEED_GENES.join(","))}&` +
    `cellTypes=${encodeURIComponent(WMG_SEED_CELL_TYPES.join(","))}&` +
    "ver=2",
  genes: WMG_SEED_GENES,
  cellTypes: WMG_SEED_CELL_TYPES,
  cellTypeIds: WMG_SEED_CELL_TYPE_IDS,
};

const WMG_SEED_TISSUES = ["blood", "lung"];
const WMG_SEED_TISSUE_IDS = ["UBERON:0000178", "UBERON:0002048"];

export const WMG_WITH_SEEDED_GENES_AND_TISSUES = {
  URL:
    `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}?` +
    `genes=${encodeURIComponent(WMG_SEED_GENES.join(","))}&` +
    `tissues=${encodeURIComponent(WMG_SEED_TISSUE_IDS.join(","))}&` +
    "ver=2",
  genes: WMG_SEED_GENES,
  tissues: WMG_SEED_TISSUES,
  tissueIds: WMG_SEED_TISSUE_IDS,
};

const WAIT_FOR_RESPONSE_TIMEOUT_MS = 10 * 1000;

/**
 * (thuang): `page.waitForResponse` sometimes times out, so we need to retry
 */
export async function goToWMG(page: Page, url?: string) {
  const targetUrl = url || `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`;
  return await tryUntil(
    async () => {
      await Promise.all([
        page.waitForResponse(
          (response) => {
            if (response.url().includes("primary_filter_dimensions")) {
              if (response.ok()) {
                return true;
              } else {
                throw new Error(
                  `Response status code is ${response.status()} for ${response.url()}`
                );
              }
            }

            return false;
          },
          { timeout: WAIT_FOR_RESPONSE_TIMEOUT_MS }
        ),
        page.goto(targetUrl),
      ]);

      await tryUntil(
        async () => {
          const numberOfTissuesBefore = await countLocator(
            page.getByTestId(TISSUE_NAME_LABEL_CLASS_NAME)
          );
          expect(numberOfTissuesBefore).toBeGreaterThan(0);
        },
        { page, silent: true }
      );
    },
    { page }
  );
}

export async function goToWMGWithSeededState(page: Page) {
  await goToWMG(page, WMG_WITH_SEEDED_GENES.URL);
  await waitForHeatmapToRender(page);
  await expandTissue(page, "lung");
  await expandTissue(page, "blood");
}

export async function waitForHeatmapToRender(page: Page) {
  await tryUntil(
    async () => {
      await expect(page.locator("canvas")).not.toHaveCount(0);
    },
    { page }
  );
}

export async function searchAndAddTissue(page: Page, tissueName: string) {
  await tryUntil(
    async () => {
      await page.getByTestId(ADD_TISSUE_ID).getByRole("button").first().click();
      await page.getByRole("tooltip").waitFor();
      await page.getByRole("tooltip").getByText(tissueName).first().click();
      // close dropdown
      await page.keyboard.press("Escape");
      await expect(
        page.getByTestId(ADD_TISSUE_ID).getByText(tissueName)
      ).toBeVisible();
    },
    { page }
  );
}

export async function addTissuesAndGenes(
  page: Page,
  tissueNames: string[],
  genes: string[]
) {
  for (const tissueName of tissueNames) {
    await searchAndAddTissue(page, tissueName);
  }
  for await (const gene of genes) {
    await Promise.all([
      // wait till loading is complete
      waitForLoadingSpinnerToResolve(page),
      searchAndAddGene(page, gene),
    ]);
  }
}
export const selectSecondaryFilterOption = async (
  page: Page,
  filterName: string
) => {
  await deleteChips({ page, filterName });

  await page.getByTestId(filterName).getByRole("button").first().click();

  // select the first option
  await selectFirstOption(page);

  // wait till loading is complete
  await page.getByText("Loading").first().waitFor({ state: "hidden" });

  // close the secondary filter pop-up
  await page.keyboard.press("Escape");

  const filterLabel = await page
    .getByTestId(filterName)
    .locator(".MuiChip-label")
    .first();

  // expect the selected filter chip to be visible under the dropdown button
  await expect(filterLabel).toBeVisible();
};

export const pickOptions = async (page: Page, n: number) => {
  for (let i = 0; i < n; i++) {
    // select the nth option
    await page.locator(`[data-option-index="${i}"]`).click();
  }
};

export const deSelectSecondaryFilterOption = async (
  page: Page,
  filterName: string
) => {
  const chipLocator = page
    .getByTestId(filterName)
    .locator(".MuiChip-deletable");

  const filterChip = await chipLocator.first();

  await expect(filterChip).toBeVisible();

  await deleteChips({ page, filterName });
};

export async function deleteChips({
  page,
  filterName,
}: {
  page: Page;
  filterName: string;
}) {
  const deleteChipLocator = page
    .getByTestId(filterName)
    .getByTestId("ClearIcon");

  await tryUntil(
    async () => {
      const deleteChips = await deleteChipLocator.all();

      deleteChips.forEach(async (deleteChip) => {
        await tryUntil(
          async () => {
            // click the delete button
            await deleteChip.click();
          },
          { page }
        );
      });

      expect(await deleteChipLocator.count()).toBe(0);
    },
    { page }
  );
}

export const selectTissueAndGeneOption = async (page: Page) => {
  // click Tissue button
  await page.getByTestId(ADD_TISSUE_ID).click();

  //pick the first 2 elements in tissue
  await pickOptions(page, 2);

  // close the pop-up
  await page.keyboard.press("Escape");

  //wait for heatmap to be visible the click action
  await page.locator('[id="heatmap-container-id"]').waitFor();

  // click Gene button
  await page.getByTestId(ADD_GENE_BTN).click();

  //pick the first n elements in tissue
  await pickOptions(page, 3);

  // close the pop-up
  await page.keyboard.press("Escape");

  //wait for gene label to appear
  await page.locator("[data-testid='gene-label-TSPAN6']").waitFor();

  //wait till loading is complete
  await page.locator(getText("Loading")).waitFor({ state: "hidden" });
};

export const checkSourceData = async (page: Page) => {
  const sourceDataButton = page.getByTestId(SOURCE_DATA_BUTTON_ID);
  const sourceDataList = page.getByTestId(SOURCE_DATA_LIST_ID);

  // click on source data icon
  await sourceDataButton.click();

  // number of element displayed on source data
  const n = await sourceDataList.locator("a").count();

  // close the pop-up
  await tryUntil(
    async () => {
      if (await sourceDataList.isVisible()) {
        await sourceDataButton.click({ force: true });
        await sourceDataList.waitFor({ state: "hidden" });
      }

      // Final check after trying to close
      const stillVisible = await sourceDataList.isVisible();
      expect(stillVisible).toBeFalsy();
    },
    {
      page,
      timeoutMs: 2 * 1000,
    }
  );

  return n;
};

export const checkPlotSize = async (page: Page) => {
  const dotLocator = page.locator('[data-zr-dom-id*="zr"]');

  await tryUntil(
    async () => {
      expect(await dotLocator.first().getAttribute("height")).toBeTruthy();
    },
    { page }
  );

  //get the number of rows on the data plot
  return page.locator('[data-zr-dom-id*="zr"]').count();
};

export async function searchAndAddGene(page: Page, geneName: string) {
  await page.getByPlaceholder(ADD_GENE_SEARCH_PLACEHOLDER_TEXT).type(geneName);
  await page.getByText(geneName).click();

  // close dropdown
  await page.keyboard.press("Escape");
}

export async function searchAndAddFilterCellType(page: Page, cellType: string) {
  const beforeCellTypeNames = await getCellTypeNames(page);
  await page.getByPlaceholder(CELL_TYPE_SEARCH_PLACEHOLDER_TEXT).fill(cellType);
  await page.getByText(cellType, { exact: true }).click();
  await page.keyboard.press("Escape");
  const afterCellTypeNames = await getCellTypeNames(page);
  expect(afterCellTypeNames.length).toBeGreaterThan(beforeCellTypeNames.length);
}
export async function removeFilteredCellType(page: Page, cellType: string) {
  const beforeCellTypeNames = await getCellTypeNames(page);
  const cellTypeTag = page.getByTestId(`cell-type-tag-${cellType}`);
  const deleteIcon = cellTypeTag.locator("svg");
  await deleteIcon.click();
  const afterCellTypeNames = await getCellTypeNames(page);
  expect(afterCellTypeNames.length).toBeLessThan(beforeCellTypeNames.length);
}
1;
