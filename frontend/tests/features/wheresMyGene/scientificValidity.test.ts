import { expect, test } from "@playwright/test";
import { getById } from "tests/utils/selectors";
import {
  goToWMG,
  searchAndAddGene,
  searchAndAddTissue,
} from "tests/utils/wmgUtil";
const { describe } = test;

describe("Rankit value tests", () => {
  test.only("Should verify MALAT1 is expressed in all values", async ({
    page,
  }) => {
    const GENE = "MALAT1";
    const TISSUE = "blood";
    await goToWMG(page);
    // add tissue and gene
    await searchAndAddTissue(page, TISSUE);
    await searchAndAddGene(page, GENE);

    // take screenshot and compare
    const screenshotTarget = page.locator(getById("blood-chart"));
    await expect(screenshotTarget).toHaveScreenshot("bloodChart.png");
  });

  test.only("Should verify is in range 0 - 7", async ({
    page,
  }) => {
    const GENE = "MALAT1";
    const TISSUE = "blood";
    await goToWMG(page);
    // add tissue and gene
    await searchAndAddTissue(page, TISSUE);
    await searchAndAddGene(page, GENE);
    await page.locator(getById("blood-chart")).hover();
    await expect(page.getByText("Tissue Composition")).toBeVisible();
  });
});