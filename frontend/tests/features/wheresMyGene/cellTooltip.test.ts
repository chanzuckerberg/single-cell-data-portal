import { expect, test } from "@playwright/test";
import { conditionallyRunTests, goToWMG } from "../../utils/wmgUtils";

import { SHARED_LINK_FILTER } from "tests/common/constants";
import { getTestID } from "tests/utils/selectors";
const { describe } = test;
describe("cell tooltip", () => {
  conditionallyRunTests();

  test(`Should verify cell tooltip hover`, async ({ page }) => {
    // set app state
    const cellName = "cell-type-name";
    await goToWMG(page, SHARED_LINK_FILTER);

    // wait for cell names to load
    await page.getByTestId(cellName).nth(0).waitFor();
    // get the number of cell names displayed
    const cells = await page.getByTestId(cellName).allInnerTexts();

    let checkedCells = 0;

    //no need to check all truncated cell names
    for (let i = 0; i < cells.length && checkedCells < 2; i++) {
      // skip cells that are not truncated
      if (!cells[i].includes("...")) continue;

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
        page.locator(getTestID("cell-type-name-tooltip")).getByText(hoverText)
      ).toBeVisible();

      checkedCells++;
    }
  });
});
