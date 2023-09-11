import { Page, expect } from "@playwright/test";
import { goToWMG } from "./wmgUtils";

const FMG_EXCLUDE_TISSUES = ["blood"];
const CELL_COUNT_ID = "cell-count";
const CELL_TYPE_NAME_ID = "cell-type-name";
const MARKER_GENE_BUTTON_ID = "marker-gene-button";
const REGEX = /^\d+\.?\d{0,2}$/;

export const ADD_GENE_SEARCH_PLACEHOLDER_TEXT = "Add Genes";

export async function addGene(page: Page, geneName: string) {
  await page.getByPlaceholder(ADD_GENE_SEARCH_PLACEHOLDER_TEXT).type(geneName);
  await page.getByText(geneName).click();

  // close dropdown
  await page.keyboard.press("Escape");
}

export async function searchAndAddGene(page: Page, geneName: string) {
  await goToWMG(page);
  await addGene(page, geneName);
}
export async function verifyAddedTissue(page: Page, tissue: string) {
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
    expect(COUNT?.replace(/\D/g, "")).toMatch(REGEX);
  }
}

export async function verifyAddedGene(page: Page, geneName: string) {
  // selected gene should be visible
  expect(await page.getByTestId(`gene-name-${geneName}`).textContent()).toBe(
    geneName
  );

  await page.getByTestId(`gene-name-${geneName}`).hover();

  // info icon
  await expect(page.getByTestId(`gene-info-icon-${geneName}`)).toBeVisible();

  // delete button
  await expect(page.getByTestId(`gene-delete-icon-${geneName}`)).toBeVisible();
}
