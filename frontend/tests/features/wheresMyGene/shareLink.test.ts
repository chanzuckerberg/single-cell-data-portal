import { Page, expect, test } from "@playwright/test";
import { LATEST_SHARE_LINK_VERSION } from "src/views/WheresMyGene/components/GeneSearchBar/components/ShareButton/utils";
import { TEST_URL } from "tests/common/constants";
import { addTissuesAndGenes, goToWMG } from "tests/utils/geneUtils";

const { describe } = test;

describe("Share link tests", () => {
  test("Should share link with single tissue and single gene", async ({
    page,
    browserName,
  }) => {
    skipFirefox(browserName);

    await setupStateAndVerifyShareLink(page, ["blood"], ["SCYL3"]);
  });

  test("Should share link with multiple tissues and multiple genes", async ({
    page,
    browserName,
  }) => {
    skipFirefox(browserName);

    await setupStateAndVerifyShareLink(
      page,
      ["blood", "lung", "liver"],
      ["SCYL3", "TSPAN6", "TNMD"]
    );
  });

  test("Should rendered shared link", async ({ page, browserName }) => {
    skipFirefox(browserName);

    const tissues = ["blood", "lung", "liver"];
    const genes = ["SCYL3", "TSPAN6", "TNMD"];

    await setupStateAndCopyShareLink(page, tissues, genes);

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

async function setupStateAndCopyShareLink(
  page: Page,
  tissues: string[],
  genes: string[]
) {
  await goToWMG(page);
  await expect(page.getByTestId("share-button")).toBeDisabled();

  // add tissues and genes
  await addTissuesAndGenes(page, tissues, genes);
  // copy share link
  await page.getByTestId("share-button").click();
}

async function setupStateAndVerifyShareLink(
  page: Page,
  tissues: string[],
  genes: string[]
) {
  await setupStateAndCopyShareLink(page, tissues, genes);

  // verify link
  await verifyShareLink(page, tissues, genes);
}

async function verifyShareLink(page: Page, tissues: string[], genes: string[]) {
  // copy link to clipboard
  const clipboardText: string = await page.evaluate(
    "navigator.clipboard.readText()"
  );

  // verify link contains all tissues, genes and version number
  const urlParams = new URLSearchParams(
    // (thuang): We only want the query params part of the URL, so we split by "?"
    decodeURIComponent(clipboardText.split("?")[1])
  );

  expect(tissues).toMatchObject(urlParams.get("tissues")?.split(",") || {});
  expect(genes).toMatchObject(urlParams.get("genes")?.split(",") || {});

  expect(urlParams.get("ver")).toBe(LATEST_SHARE_LINK_VERSION);

  // verify encoding
  // encode tissues and genes first
  const encodedTissues = encodeURIComponent(tissues.toString());
  const encodedGenes = encodeURIComponent(genes.toString());

  const link = `${TEST_URL}/gene-expression?tissues=${encodedTissues}&genes=${encodedGenes}&ver=${LATEST_SHARE_LINK_VERSION}`;

  expect(clipboardText).toBe(link);
}

function skipFirefox(browserName: string) {
  test.skip(browserName === "firefox", "No Clipboard read permission");
}
