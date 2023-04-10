import { test } from "@playwright/test";
import { ADD_TISSUE_BTN } from "tests/utils/constants";
import { goToWMG, verifyAddedTissue } from "tests/utils/geneUtils";

const { describe } = test;

describe("Add tissue tests", () => {
  test("Should select tissue using keyboard arrow key to select", async ({
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

  test("Should select tissue by searching", async ({ page }) => {
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
