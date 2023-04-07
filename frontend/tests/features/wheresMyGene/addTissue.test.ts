import { expect, Page, test } from "@playwright/test";
import { TEST_URL } from "../../common/constants";
import { ROUTES } from "src/common/constants/routes";
import {
  ADD_GENE_BTN,
  ADD_GENE_LBL,
  ADD_TISSUE_BTN,
  ADD_TISSUE_LBL,
} from "tests/utils/constants";
const REGEX = /^\d+\.?\d{0,2}$/;
const { describe } = test;

function goToWMG(page: Page) {
  return Promise.all([
    page.waitForResponse(
      (resp: { url: () => string | string[]; status: () => number }) =>
        resp.url().includes("/wmg/v1/filters") && resp.status() === 200
    ),
    page.goto(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`),
  ]);
}
describe("Add tissue tests", () => {
  test.only("Should select tissue by keyboard down key", async ({ page }) => {
    await goToWMG(page);
    // click +Tissue button
    await page.getByTestId(ADD_TISSUE_BTN).click();
    await page.keyboard.press("ArrowDown");
    await page.keyboard.press("ArrowDown");

    // select 2nd element
    await page.keyboard.press("Enter");

    // close dropdown
    await page.keyboard.press("Escape");
    await expect(page.getByTestId("cell-type-labels-lung")).toBeVisible();
    const CELL_COUNTS = page.getByTestId("cell-type-label-count");
    expect(await CELL_COUNTS.count()).toBe(151); // flaky, will break with changes in test data
    for (let i = 0; i < (await CELL_COUNTS.count()); i++) {
      // cell name
      expect(
        CELL_COUNTS.nth(i).getByTestId("cell-type-name").textContent()
      ).not.toBeUndefined();

      // info icon
      console.log(
        await CELL_COUNTS.nth(i).getByTestId("cell-type-name").textContent()
      );
      expect(
        CELL_COUNTS.nth(i).getByTestId("marker-gene-button")
      ).toBeVisible();
      // cell count
      const COUNT = await CELL_COUNTS.nth(i)
        .getByTestId("cell-count")
        .textContent();
      expect(COUNT?.replace("k", "")).toMatch(REGEX);
    }
  });
});
