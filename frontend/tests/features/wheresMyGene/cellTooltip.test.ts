import { expect, test } from "@playwright/test";
import {
  conditionallyRunTests,
  goToWMGWithSeededState,
} from "../../utils/wmgUtils";
import { tryUntil } from "tests/utils/helpers";
const { describe } = test;

describe("cell tooltip", () => {
  /**
   * TODO(thuang): Remove forceRun when all WMG e2e tests are enabled.
   * `forceRun` is just to incrementally add tests back in the meantime
   */
  conditionallyRunTests({ forceRun: true });

  test(`Should verify cell tooltip hover`, async ({ page }) => {
    // set app state
    const cellName = "cell-type-name";
    await goToWMGWithSeededState(page);

    // wait for cell names to load
    await page.getByTestId(cellName).nth(0).waitFor();

    // get the number of cell names displayed
    const beforeCells = await page.getByTestId(cellName).allInnerTexts();

    // (thuang): Expand blood tissue to find truncated cell names
    await tryUntil(
      async () => {
        const afterCells = await page.getByTestId(cellName).allInnerTexts();
        expect(afterCells.length).toBeGreaterThan(beforeCells.length);
      },
      { page }
    );

    // get the number of cell names displayed
    const cells = await page.getByTestId(cellName).allInnerTexts();

    let checkedCells = 0;

    // no need to check all truncated cell names
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

      await tryUntil(
        async () => {
          await page.getByTestId(cellName).nth(i).hover();

          //expect the hover text to be displayed
          await expect(page.getByText(hoverText)).toBeVisible();
        },
        { page }
      );

      checkedCells++;
    }
  });
});
