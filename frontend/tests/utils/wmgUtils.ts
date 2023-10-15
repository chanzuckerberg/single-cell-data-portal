import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "../common/constants";
import { expect, Page } from "@playwright/test";
import { getTestID, getText } from "tests/utils/selectors";
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

  const filter_label = `${getTestID(filterName)} [role="button"]`;

  // expect the selected filter chip to be visible under the dropdown button
  await expect(page.locator(filter_label).first()).toBeVisible();
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
  //click on source data icon
  await page.locator('[data-testid="source-data-button"]').click();

  // number of element displayed on source data
  const n = await page.locator('[data-testid="source-data-list"] a').count();

  // close the pop-up
  /**
   * (thuang): Sometimes pressing escape once wasn't closing the side panel, so
   * wrapping this to retry and assert the panel is indeed closed
   */
  await tryUntil(
    async () => {
      await page.keyboard.press("Escape");

      await tryUntil(
        async () => {
          expect(
            await page.locator('[data-testid="source-data-list"]').isVisible()
          ).toBeFalsy();
        },
        {
          page,
          /**
           * (thuang): we don't need to wait for too long to retry pressing escape
           * button, since the source data panel should close within 2s
           */
          maxRetry: 10,
        }
      );
    },
    { page }
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
  const deleteIcon = cellTypeTag.getByTestId("CancelIcon");
  await deleteIcon.click();
  const afterCellTypeNames = await getCellTypeNames(page);
  expect(afterCellTypeNames.length).toBeLessThan(beforeCellTypeNames.length);
}
1;
