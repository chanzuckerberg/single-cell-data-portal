import { expect, Page, test } from "@playwright/test";
import { TEST_URL } from "../../common/constants";
import { ROUTES } from "src/common/constants/routes";
import { ADD_TISSUE_BTN, ADD_TISSUE_LBL } from "tests/utils/constants";
const REGEX = /^\d+\.?\d{0,2}$/;
const { describe } = test;
const FMG_EXCLUDE_TISSUES = ["blood"];
const CELL_COUNT_ID = "cell-count";
const CELL_TYPE_NAME_ID = "cell-type-name";
const MARKER_GENE_BUTTON_ID = "marker-gene-button";

function goToWMG(page: Page) {
  return Promise.all([
    page.waitForResponse(
      (resp: { url: () => string | string[]; status: () => number }) =>
        resp.url().includes("/wmg/v1/filters") && resp.status() === 200
    ),
    page.goto(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`),
  ]);
}
async function verifyAddedTissue(page: Page, tissue: string) {
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

describe("Add tissue tests", () => {
  test.only("Should select tissue using keyboard arrow key to select", async ({
    page,
  }) => {
    const TISSUE = "lung";
    await goToWMG(page);
    // click +Tissue button
    await page.getByTestId(ADD_TISSUE_BTN).click();
    await page.keyboard.press("ArrowDown");
    await page.keyboard.press("ArrowDown");

    // select 2nd element
    await page.keyboard.press("Enter");

    // close dropdown
    await page.keyboard.press("Escape");

    // verify selected tissue details
    await verifyAddedTissue(page, TISSUE);
  });

  test.only("Should select tissue by searching", async ({ page }) => {
    const TISSUE = "blood";
    await goToWMG(page);
    // click +Tissue button
    await page.getByTestId(ADD_TISSUE_BTN).click();
    await page.getByPlaceholder("Search").type(TISSUE);
    await page.getByText(TISSUE).click();

    // close dropdown
    await page.keyboard.press("Escape");
    // verify selected tissue details
    await verifyAddedTissue(page, TISSUE);
  });
});
