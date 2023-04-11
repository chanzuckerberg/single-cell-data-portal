import { expect, test } from "@playwright/test";
import { ADD_GENE_BTN } from "tests/utils/constants";
import { goToWMG, verifyAddedGene } from "tests/utils/geneUtils";
import uaParser from "ua-parser-js";
const { describe } = test;

describe("Add gene tests", () => {
  test("Should select tissue using keyboard arrow key to select", async ({
    page,
  }) => {
    const GENE = "TNMD";
    await goToWMG(page);
    // click +Tissue button
    await page.getByTestId(ADD_GENE_BTN).click();
    await page.keyboard.press("ArrowDown");
    await page.keyboard.press("ArrowDown");

    // select 2nd element
    await page.keyboard.press("Enter");

    // close dropdown
    await page.keyboard.press("Escape");

    // verify selected tissue details
    await verifyAddedGene(page, GENE);
  });

  test("Should select gene by searching", async ({ page }) => {
    const GENE = "FGR";
    await goToWMG(page);
    // click +Tissue button
    await page.getByTestId(ADD_GENE_BTN).click();
    await page.getByPlaceholder("Search").type(GENE);
    await page.getByText(GENE).click();

    // close dropdown
    await page.keyboard.press("Escape");
    // verify selected tissue details
    await verifyAddedGene(page, GENE);
  });

  test("Should select gene by copy pasting", async ({ page, browserName }) => {
    test.skip(browserName === "firefox", "No Clipboard read permission");
    const TEST_GENES = "DMP1,SCYL3,CFH";
    await goToWMG(page);
    // click +Tissue button
    await page.getByTestId(ADD_GENE_BTN).click();

    // copy & paste clipboard contents into search
    // we will first write into search field so we have what to copy & paste

    await page.getByPlaceholder("Search").type(TEST_GENES);
    const getUA = await page.evaluate(() => navigator.userAgent);
    const userAgentInfo = uaParser(getUA);
    const modifier = userAgentInfo.os.name?.includes("Mac")
      ? "Meta"
      : "Control";
    await page.getByPlaceholder("Search").focus();
    await page.keyboard.press(`${modifier}+KeyA`);
    await page.keyboard.press(`${modifier}+KeyC`);
    await page.getByPlaceholder("Search").click();
    await page.keyboard.press(`${modifier}+KeyV`);
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
    await goToWMG(page);
    // click +Tissue button
    await page.getByTestId(ADD_GENE_BTN).click();
    await page.getByPlaceholder("Search").type(GENE);
    await page.getByText(GENE).click();

    // close dropdown
    await page.keyboard.press("Escape");

    // delete gene
    await page.getByTestId(`gene-name-${GENE}`).hover();
    await page.getByTestId(`gene-delete-icon-${GENE}`).click();
    // check gene is deleted
    expect(page.getByTestId(`gene-name-${GENE}`)).not.toBeVisible();
  });
});
