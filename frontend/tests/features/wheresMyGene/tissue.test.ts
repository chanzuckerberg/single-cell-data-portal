import { test } from "@playwright/test";
import { ADD_TISSUE_BTN } from "tests/utils/constants";
import {
  goToWMG,
  searchAndAddTissue,
  verifyAddedTissue,
} from "tests/utils/geneAndTissueUtils";

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
    await searchAndAddTissue(page, TISSUE);
    // verify selected tissue details
    await verifyAddedTissue(page, TISSUE);
  });
});
