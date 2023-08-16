import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "../common/constants";
import { expect, Page, test } from "@playwright/test";
import { getTestID, getText } from "tests/utils/selectors";
import { expandTissue, selectFirstOption, tryUntil } from "./helpers";
import { ADD_GENE_BTN, ADD_TISSUE_ID } from "../common/constants";
import { ADD_GENE_SEARCH_PLACEHOLDER_TEXT } from "tests/utils/geneUtils";

const { skip, beforeEach } = test;

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

const ENVS_TO_RUN_TESTS = [
  "api.cellxgene.dev.single-cell.czi.technology",
  "api.cellxgene.staging.single-cell.czi.technology",
  "api.cellxgene.cziscience.com",
];

export function conditionallyRunTests({
  forceRun = false,
}: {
  forceRun?: boolean;
} = {}) {
  if (forceRun) return;

  /**
   * (thuang): `beforeEach()` is needed, since without it, `describe()` will
   * just eager inventory tests to run BEFORE our global setup sets `process.env.API_URL`
   */
  beforeEach(() => {
    skip(
      // (thuang): Temporarily skip WMG tests
      // (thuang): Temporarily skip WMG tests
      // (thuang): Temporarily skip WMG tests
      ENVS_TO_RUN_TESTS.every(() => true),
      // ENVS_TO_RUN_TESTS.every((env) => !process.env.API_URL?.includes(env)),
      "WMG tests only work with dev/staging/prod API URLs"
    );
  });
}

/**
 * (thuang): `page.waitForResponse` sometimes times out, so we need to retry
 */
export async function goToWMG(page: Page, url?: string) {
  const targetUrl = url || `${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`;
  return await tryUntil(
    async () => {
      await Promise.all([
        page.waitForResponse(
          (resp: { url: () => string | string[]; status: () => number }) =>
            resp.url().includes("/wmg/v2/filters") && resp.status() === 200
        ),
        page.goto(targetUrl),
      ]);
    },
    { page }
  );
}

export async function goToWMGWithSeededState(page: Page) {
  await goToWMG(page, WMG_WITH_SEEDED_GENES.URL);
  await expandTissue(page, "lung");
  await expandTissue(page, "blood");
  await waitForHeatmapToRender(page);
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
    },
    { page }
  );

  await page.getByRole("tooltip").getByText(tissueName).first().click();

  // close dropdown
  await page.keyboard.press("Escape");
}

export async function addTissuesAndGenes(
  page: Page,
  tissueNames: string[],
  genes: string[]
) {
  for await (const tissueName of tissueNames) {
    await searchAndAddTissue(page, tissueName);
  }
  for await (const gene of genes) {
    await Promise.all([
      // wait till loading is complete
      page.getByText("Loading").first().waitFor({ state: "hidden" }),
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
