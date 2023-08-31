import { expect, test } from "@playwright/test";
import {
  ADD_GENE_SEARCH_PLACEHOLDER_TEXT,
  goToWMG,
  searchAndAddGene,
  verifyAddedGene,
} from "tests/utils/geneUtils";
import { selectNthOption } from "tests/utils/helpers";
import { conditionallyRunTests } from "tests/utils/wmgUtils";
const { describe } = test;

describe("Manage gene tests", () => {
  /**
   * TODO(thuang): Remove forceRun when all WMG e2e tests are enabled.
   * `forceRun` is just to incrementally add tests back in the meantime
   */
  conditionallyRunTests({ forceRun: true });

  test("Should select gene using keyboard arrow key to select", async ({
    page,
  }) => {
    const GENE = "TNMD";
    await goToWMG(page);
    // click +Gene button
    await page.getByPlaceholder(ADD_GENE_SEARCH_PLACEHOLDER_TEXT).click();

    // use arrow down buttons to select the 2nd option
    await selectNthOption(page, 3);

    // verify selected gene details
    await verifyAddedGene(page, GENE);
  });

  test("Should select gene by searching", async ({ page }) => {
    const GENE = "FGR";
    await searchAndAddGene(page, GENE);

    // verify selected tissue details
    await verifyAddedGene(page, GENE);
  });

  test("Should select gene by comma separated list", async ({ page }) => {
    const TEST_GENES = "DMP1,SCYL3,CFH";

    await goToWMG(page);

    await page
      .getByPlaceholder(ADD_GENE_SEARCH_PLACEHOLDER_TEXT)
      .type(TEST_GENES);

    await page.keyboard.press("Enter");

    // verify selected tissue details
    // 'for' loop helps avoid race condition
    const GENES = TEST_GENES.split(",");
    for (let i = 0; i < GENES.length; i++) {
      await verifyAddedGene(page, GENES[i]);
    }
  });

  test("Should remove gene", async ({ page }) => {
    const GENE = "SCYL3";
    await searchAndAddGene(page, GENE);

    // delete gene
    await page.getByTestId(`gene-name-${GENE}`).hover();
    await page.getByTestId(`gene-delete-icon-${GENE}`).click();

    // check gene is deleted
    expect(page.getByTestId(`gene-name-${GENE}`)).not.toBeVisible();
  });
});
