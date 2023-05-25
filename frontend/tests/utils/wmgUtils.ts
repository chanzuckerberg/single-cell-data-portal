import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "../common/constants";
import { expect, Page } from "@playwright/test";
import { getTestID, getText } from "tests/utils/selectors";
import { selectFirstOption, tryUntil } from "./helpers";
import { ADD_GENE_BTN, ADD_TISSUE_BTN } from "../common/constants";

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
export const selectSecondaryFilterOption = async (
  page: Page,
  filterName: string
) => {
  await page.getByTestId(filterName).getByRole("button").click();

  // select the first option
  await selectFirstOption(page);

  // close the secondary filter pop-up
  await page.keyboard.press("Escape");

  const filter_label = `${getTestID(filterName)} [role="button"]`;
  // expect the selected filter chip to be visible under the dropdown button
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

export const deSelectSecondaryFilterOption = async (
  page: Page,
  filterName: string
) => {
  const filter_label = `${getTestID(filterName)} [role="button"]`;
  // expect the selected filter to be visible
  await expect(page.locator(filter_label)).toBeVisible();

  // click the cancel button
  await page.getByTestId("ClearIcon").click();

  // verify the selected filter is not visible
  const visibility = await page.locator(filter_label).isVisible();
  expect(visibility).toBeFalsy();
};

export const selectTissueAndGeneOption = async (page: Page) => {
  // click Tissue button
  await page.getByTestId(ADD_TISSUE_BTN).click();

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

export async function searchAndAddGene(page: Page, geneName: string) {
  // click +Tissue button
  await page.getByTestId(ADD_GENE_BTN).click();
  await page.getByPlaceholder("Search").type(geneName);
  await page.getByText(geneName).click();

  // close dropdown
  await page.keyboard.press("Escape");
}
