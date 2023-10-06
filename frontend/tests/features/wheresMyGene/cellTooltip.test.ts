import { expect, test } from "@playwright/test";
import { expandTissue, getCellTypeNames, goToPage } from "tests/utils/helpers";
import { TEST_URL } from "tests/common/constants";
import { ROUTES } from "src/common/constants/routes";
import {
  CELL_TYPE_NAME_LABEL_CLASS_NAME,
  CELL_TYPE_ROW_CLASS_NAME,
} from "src/views/WheresMyGeneV2/components/HeatMap/components/YAxisChart/constants";
const { describe } = test;

describe("cell tooltip", () => {
  test(`Should verify cell tooltip hover`, async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`, page);

    const cellName = CELL_TYPE_NAME_LABEL_CLASS_NAME;

    // (thuang): Expand blood tissue to find truncated cell names
    await expandTissue(page, "blood");

    const cellsAfter = await getCellTypeNames(page);

    // get the number of cell names displayed

    let checkedCells = 0;

    // no need to check all truncated cell names
    for (let i = 0; i < cellsAfter.length && checkedCells < 2; i++) {
      // skip cells that are not truncated
      if (!cellsAfter[i].includes("...")) continue;

      // get the text displayed when user hovers over truncated cell names
      const hoverText =
        (await page
          .getByTestId(CELL_TYPE_ROW_CLASS_NAME)
          .nth(i)
          .getByTestId("cell-type-full-name")
          .textContent()) || "no text";

      await page.getByTestId(cellName).nth(i).hover();

      //expect the hover text to be displayed
      await expect(page.getByText(hoverText)).toBeVisible();

      checkedCells++;
    }
  });
});
