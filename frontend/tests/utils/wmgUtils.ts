import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "../common/constants";
import { expect, Page } from "@playwright/test";
import { getTestID, getText } from "tests/utils/selectors";
import { tryUntil } from "./helpers";
import {
  ADD_GENE_BTN,
  ADD_TISSUE_BTN,
  TWO_DECIMAL_NUMBER_REGEX,
} from "../common/constants";

const FMG_EXCLUDE_TISSUES = ["blood"];
const CELL_COUNT_ID = "cell-count";
const CELL_TYPE_NAME_ID = "cell-type-name";
const MARKER_GENE_BUTTON_ID = "marker-gene-button";

export const selectFilterOption = async (page: Page, filterName: string) => {
  // click the filter at the corner this is done due to the fact that the default click is being intercepted by another element
  await page.getByTestId(filterName).getByRole("button").click();

  // select the first option
  await page.locator("[data-option-index='0']").click();

  // close the pop-up
  await page.getByTestId("dataset-filter").click();

  const filter_label = `${getTestID(filterName)} [role="button"]`;
  // expect the selected filter to be visible
  await expect(page.locator(filter_label)).toBeVisible();

  //wait till loading is complete
  await page.locator(getText("Loading")).waitFor({ state: "hidden" });
};
export const pickOptions = async (page: Page, n: number) => {
  for (let i = 0; i < n; i++) {
    // select the nth option
    await page.locator(`[data-option-index="${i}"]`).click();
  }
};

export const deSelectFilterOption = async (page: Page, filterName: string) => {
  const filter_label = `${getTestID(filterName)} [role="button"]`;
  // expect the selected filter to be visible
  await expect(page.locator(filter_label)).toBeVisible();

  // click the cancel button
  await page.getByTestId("ClearIcon").click();

  // verify the selected filter is not visible
  const visibility = await page.locator(filter_label).isVisible();
  expect(visibility).toBeFalsy();
};

export const selectOption = async (page: Page, filterName: string) => {
  // click the filter
  await page.getByTestId(filterName).click();
};

export const selectTissueAndGeneOption = async (page: Page) => {
  // click Tissue button
  await selectOption(page, "add-tissue-btn");

  //pick the first 2 elements in tissue
  await pickOptions(page, 2);

  // close the pop-up
  await page.keyboard.press("Escape");

  //wait for heatmap to be visible the click action
  await page.locator('[id="heatmap-container-id"]').waitFor();

  // click Gene button
  await selectOption(page, "add-gene-btn");

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
  //get the number of rows on the data plot
  const n = await page.locator('[data-zr-dom-id*="zr"]').count();
  let sumOfHeights = 0;
  for (let i = 0; i < n; i++) {
    const row = await page.locator('[data-zr-dom-id*="zr"]').nth(i);

    const height = await row.getAttribute("height");

    if (height !== null) {
      sumOfHeights += parseInt(height);
    }
  }
  return sumOfHeights;
};

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
            resp.url().includes("/wmg/v1/filters") && resp.status() === 200
        ),
        page.goto(targetUrl),
      ]);
    },
    { page }
  );
}
export async function searchAndAddGene(page: Page, geneName: string) {
  // click +Tissue button
  await page.getByTestId(ADD_GENE_BTN).click();
  await page.getByPlaceholder("Search").type(geneName);
  await page.getByText(geneName).click();

  // close dropdown
  await page.keyboard.press("Escape");
}
export async function searchAndAddTissue(page: Page, tissueName: string) {
  // click +Tissue button
  await page.getByTestId(ADD_TISSUE_BTN).click();
  await page.getByPlaceholder("Search").type(tissueName);
  await page.getByText(tissueName).click();

  // close dropdown
  await page.keyboard.press("Escape");
}
export async function addTissuesAndGenes(
  page: Page,
  tissues: string[],
  genes: string[]
) {
  for await (const tissue of tissues) {
    await searchAndAddTissue(page, tissue);
  }
  for await (const gene of genes) {
    await searchAndAddGene(page, gene);
  }
}
export async function verifyAddedTissue(page: Page, tissue: string) {
  // STEP 1 & Add Tissues texts should disappear
  await expect(page.getByText("STEP 1")).not.toBeVisible();
  await expect(page.getByTestId("Add Tissues")).not.toBeVisible();

  // selected tissue should be visible
  await expect(page.getByTestId(`cell-type-labels-${tissue}`)).toBeVisible();

  // verify cell counts: name, icon and count
  const CELL_COUNTS = page.getByTestId("cell-type-label-count");
  for (let i = 0; i < (await CELL_COUNTS.count()); i++) {
    const COUNT = await CELL_COUNTS.nth(i)
      .getByTestId(CELL_COUNT_ID)
      .textContent();
    // cell name
    expect(
      CELL_COUNTS.nth(i).getByTestId(CELL_TYPE_NAME_ID).textContent()
    ).not.toBeUndefined();

    // info icon: if not blood and count is > 25
    if (
      !FMG_EXCLUDE_TISSUES.includes(tissue) &&
      Number(COUNT?.replace(/\D/g, "")) > 25
    ) {
      expect(
        CELL_COUNTS.nth(i).getByTestId(MARKER_GENE_BUTTON_ID)
      ).toBeVisible();
    }

    // cell count
    expect(COUNT?.replace(/\D/g, "")).toMatch(TWO_DECIMAL_NUMBER_REGEX);
  }
}

export async function verifyAddedGene(page: Page, geneName: string) {
  // STEP 1 & Add Tissues texts should disappear
  await expect(page.getByText("STEP 2")).not.toBeVisible();
  await expect(page.getByTestId("Add Genes")).not.toBeVisible();

  // selected gene should be visible
  expect(await page.getByTestId(`gene-name-${geneName}`).textContent()).toBe(
    geneName
  );

  // info icon
  await expect(page.getByTestId(`gene-info-icon-${geneName}`)).toBeVisible();

  // delete button
  await page.getByTestId(`gene-name-${geneName}`).hover();
  await expect(page.getByTestId(`gene-delete-icon-${geneName}`)).toBeVisible();
}
