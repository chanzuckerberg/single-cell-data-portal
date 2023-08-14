import { test } from "@playwright/test";
import { ADD_TISSUE_ID } from "tests/common/constants";
import { goToWMG, verifyAddedTissue } from "tests/utils/geneUtils";
import { selectNthOption, tryUntil } from "tests/utils/helpers";
import { conditionallyRunTests } from "tests/utils/wmgUtils";

const { describe } = test;

describe("Add tissue tests", () => {
  /**
   * TODO(thuang): Remove forceRun when all WMG e2e tests are enabled.
   * `forceRun` is just to incrementally add tests back in the meantime
   */
  conditionallyRunTests({ forceRun: true });

  test("Should select tissue using keyboard arrow key to select", async ({
    page,
  }) => {
    const TISSUE = "lung";
    await goToWMG(page);
    // click +Tissue button
    await page.getByTestId(ADD_TISSUE_ID).click();

    // select second option
    await selectNthOption(page, 3);

    // verify selected tissue details
    await tryUntil(
      async () => {
        await verifyAddedTissue(page, TISSUE);
      },
      { page }
    );
  });

  test("Should select tissue by searching", async ({ page }) => {
    const TISSUE = "blood";
    await goToWMG(page);
    // click +Tissue button
    await page.getByTestId(ADD_TISSUE_ID).click();
    await page.getByPlaceholder("Search").type(TISSUE);
    await page.getByText(TISSUE).click();

    // close dropdown
    await page.keyboard.press("Escape");
    // verify selected tissue details

    await tryUntil(
      async () => {
        await verifyAddedTissue(page, TISSUE);
      },
      { page }
    );
  });
});
