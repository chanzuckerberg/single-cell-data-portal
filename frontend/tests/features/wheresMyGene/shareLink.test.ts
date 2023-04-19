import { expect, test } from "@playwright/test";
import {
  goToWMG,
  searchAndAddGene,
  searchAndAddTissue,
} from "tests/utils/geneUtils";

const { describe } = test;

describe("Share link tests", () => {
  test.only("Should share link", async ({ page, browserName }) => {
    test.skip(browserName === "firefox", "No Clipboard read permission");
    const TISSUE = "blood";
    const GENE = "SCYL3";
    await goToWMG(page);
    await expect(page.getByTestId("share-button")).toBeDisabled();
    // add tissue
    await searchAndAddTissue(page, TISSUE);
    // add gene
    await searchAndAddGene(page, GENE);

    // copy share link
    await page.getByTestId("share-button").click();
    // verify content
    const clipboardText = await page.evaluate("navigator.clipboard.readText()");
    expect(clipboardText).toContain(`tissues=${TISSUE}&genes=${GENE}&ver=2`);
  });
});

//https://localhost:3000/gene-expression?tissues=blood%2CUBERON%3A0000966+%28organoid%29%2Cparacolic+gutter&genes=TNMD%2CRP5-1142A6.8%2CCAPN12&ver=2