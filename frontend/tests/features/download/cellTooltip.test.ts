import { expect, test } from "@playwright/test";
import { goToWMG } from "../../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";

import { SHARED_LINK_NO_FILTER } from "tests/common/constants";
const { describe } = test;
describe("cell tool tip", () => {
  //skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test.only(`Should verify cell tool tip hover`, async ({ page }) => {
    // set app state
    const cellName = "cell-type-name";
    await goToWMG(page, SHARED_LINK_NO_FILTER);
    await page.getByTestId(cellName).nth(0).waitFor();
    const cells = await page.getByTestId(cellName).allInnerTexts();
    for (let i = 0; i < cells.length; i++) {
      if (cells[i].includes("...")) {
        const hoverText =
          (await page
            .getByTestId("cell-type-label-count")
            .nth(i)
            .getByTestId("cell-type-full-name")
            .textContent()) || "no text";
        console.log(hoverText);
        await page.getByTestId(cellName).nth(i).hover();
        await expect(
          page.locator('[data-foo="bar"] div').getByText(hoverText)
        ).toBeVisible();
      }
    }
  });
});
