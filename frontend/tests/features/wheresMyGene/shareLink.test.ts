import { Page, expect, test } from "@playwright/test";
import { TEST_URL } from "tests/common/constants";
import { addTissuesAndGenes, goToWMG } from "tests/utils/geneUtils";

const { describe } = test;

describe("Share link tests", () => {
  test.only("Should share link with single tissue and single gene", async ({
    page,
    browserName,
  }) => {
    test.skip(browserName === "firefox", "No Clipboard read permission");
    const tissues = ["blood"];
    const genes = ["SCYL3"];
    await goToWMG(page);
    await expect(page.getByTestId("share-button")).toBeDisabled();

    // add tissues and genes
    await addTissuesAndGenes(page, tissues, genes);
    // copy share link
    await page.getByTestId("share-button").click();

    // verify link
    await verifyShareLink(page, tissues, genes);
  });
  test.only("Should share link with multiple tissues and multiple genes", async ({
    page,
    browserName,
  }) => {
    test.skip(browserName === "firefox", "No Clipboard read permission");
    const tissues = ["blood", "lung", "liver"];
    const genes = ["SCYL3", "TSPAN6", "TNMD"];
    await goToWMG(page);
    await expect(page.getByTestId("share-button")).toBeDisabled();

    // add tissues and genes
    await addTissuesAndGenes(page, tissues, genes);
    // copy share link
    await page.getByTestId("share-button").click();

    // verify link
    await verifyShareLink(page, tissues, genes);
  });
  test.only("Should rendered shared link", async ({ page, browserName }) => {
    test.skip(browserName === "firefox", "No Clipboard read permission");
    const tissues = ["blood", "lung", "liver"];
    const genes = ["SCYL3", "TSPAN6", "TNMD"];
    await goToWMG(page);
    // add tissues and genes
    await addTissuesAndGenes(page, tissues, genes);
    // copy share link
    await page.getByTestId("share-button").click();

    const clipboardText: string = await page.evaluate(
      "navigator.clipboard.readText()"
    );

    await goToWMG(page, clipboardText);
    tissues.forEach(async (tissue) => {
      // selected tissue should be visible
      await expect(
        page.getByTestId(`cell-type-labels-${tissue}`)
      ).toBeVisible();
    });
    genes.forEach(async (gene) => {
      // selected gene should be visible
      expect(await page.getByTestId(`gene-name-${gene}`).textContent()).toBe(
        gene
      );
    });
  });
});

async function verifyShareLink(page: Page, tissues: string[], genes: string[]) {
  // copy link to clipboard
  const clipboardText: string = await page.evaluate(
    "navigator.clipboard.readText()"
  );

  // verify link contains all tissues, genes and version number
  const urlParams = new URLSearchParams(clipboardText);
  urlParams.getAll("tissues").forEach((x) => expect(tissues).toContain(x));
  expect(urlParams.getAll("genes").toString()).toBe(genes.toString()); //todo: above syntax does not work for genes
  expect(urlParams.get("ver")).toBe("2");

  // verify encoding
  // encode tissues and genes first
  const encodedTissues = encodeURIComponent(tissues.toString());
  const encodedGenes = encodeURIComponent(genes.toString());
  const link = `${TEST_URL}/gene-expression?tissues=${encodedTissues}&genes=${encodedGenes}&ver=2`;
  expect(clipboardText).toBe(link);
}
