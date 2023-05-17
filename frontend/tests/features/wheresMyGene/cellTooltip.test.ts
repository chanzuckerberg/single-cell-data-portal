import { expect, test } from "@playwright/test";
import { goToWMG } from "../../utils/wmgUtils";
import { isDevStagingProd } from "tests/utils/helpers";

import { SHARED_LINK_FILTER } from "tests/common/constants";
const { describe, skip } = test;
describe("cell tool tip", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  test(`Should verify cell tool tip hover`, async ({ page }) => {
    // set app state
    const cellName = "cell-type-name";
    await goToWMG(page, SHARED_LINK_FILTER);

    // wait for cell names to load
    await page.getByTestId(cellName).nth(0).waitFor();
    // get the number of cell names displayed
    const cells = await page.getByTestId(cellName).allInnerTexts();
    for (let i = 0; i < cells.length; i++) {
      if (cells[i].includes("...")) {
        // get the text displayed when user hovers over truncated cell names
        const hoverText =
          (await page
            .getByTestId("cell-type-label-count")
            .nth(i)
            .getByTestId("cell-type-full-name")
            .textContent()) || "no text";
        await page.getByTestId(cellName).nth(i).hover();

        //expect the hover text to be displayed
        await expect(
          page.locator('[data-foo="bar"] div').getByText(hoverText)
        ).toBeVisible();
      }
      //no need to check all truncated cell names
      if (i > 2) {
        break;
      }
    }
  });
});
