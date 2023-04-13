import { Page, expect } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "tests/common/constants";
import { ADD_GENE_BTN, ADD_TISSUE_BTN } from "./constants";

const FMG_EXCLUDE_TISSUES = ["blood"];
const CELL_COUNT_ID = "cell-count";
const CELL_TYPE_NAME_ID = "cell-type-name";
const MARKER_GENE_BUTTON_ID = "marker-gene-button";
const REGEX = /^\d+\.?\d{0,2}$/;

export async function goToWMG(page: Page) {
  return Promise.all([
    page.waitForResponse(
      (resp: { url: () => string | string[]; status: () => number }) =>
        resp.url().includes("/wmg/v1/filters") && resp.status() === 200
    ),
    page.goto(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`),
  ]);
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
    expect(COUNT?.replace(/\D/g, "")).toMatch(REGEX);
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
