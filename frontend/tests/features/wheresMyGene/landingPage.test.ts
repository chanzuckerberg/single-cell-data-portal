import { expect, Page } from "@playwright/test";
import { TEST_URL } from "../../common/constants";
import { ROUTES } from "src/common/constants/routes";
import { ADD_GENE_BTN } from "tests/common/constants";
import { getById } from "tests/utils/selectors";
import { tryUntil } from "tests/utils/helpers";
import { test } from "tests/common/test";

const { describe } = test;
const ALERT =
  "Subscribeto our newsletter to receive updates about new features.Send us feedback with thisquick survey.";

const SURVEY_LINK = "https://airtable.com/shrLwepDSEX1HI6bo";

const LEGEND_WRAPPER = "legend-wrapper";
const DOT_SIZES = ["4", "9", "12", "14", "16"];

function goToWMG(page: Page) {
  return Promise.all([
    page.waitForResponse(
      (resp: { url: () => string | string[]; status: () => number }) =>
        resp.url().includes("/wmg/v2/filters") && resp.status() === 200
    ),
    page.goto(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`),
  ]);
}
describe("Tests for Gene Expression page", () => {
  test("Should verify main panel components", async ({ page }) => {
    await goToWMG(page);
    // +Gene button
    await expect(page.getByTestId(ADD_GENE_BTN)).toBeVisible();
    // survey alert
    await expect(
      page.getByTestId("newsletter-modal-banner-wrapper")
    ).toContainText(ALERT);
    await expect(page.getByText("quick survey")).toHaveAttribute(
      "href",
      SURVEY_LINK
    );
    // default organism filter
    await expect(page.getByTestId("add-organism")).toContainText(
      "Homo sapiens"
    );

    // Download
    await expect(page.getByTestId(LEGEND_WRAPPER)).toContainText("Download");
    await expect(page.getByTestId("download-button")).toBeDisabled();

    // Share
    await expect(page.getByTestId(LEGEND_WRAPPER)).toContainText("Share");
    await expect(page.getByTestId("share-button")).toBeDisabled();

    // Source data
    await expect(page.getByTestId(LEGEND_WRAPPER)).toContainText("Source Data");
    await expect(
      page.getByTestId("source-data-button").locator("svg")
    ).toBeEnabled();

    // Gene expression
    await expect(page.getByTestId(LEGEND_WRAPPER)).toContainText(
      "Gene Expression"
    );
    await expect(
      page.locator('[id="visualization-color-scale"]')
    ).toBeVisible();

    // Gene expression in cells
    await expect(page.getByTestId(LEGEND_WRAPPER)).toContainText(
      "Expressed in Cells (%)"
    );
    await expect(page.locator('[id="expressed-in-cells-dots"]')).toBeVisible();
    expect(
      await page.$$eval(
        '[data-testid="expressed-in-cells-dots-size"]',
        (dots) => dots.map((d) => d.getAttribute("size"))
      )
    ).toStrictEqual(DOT_SIZES);
  });
  test("Should verify top nav", async ({ page }) => {
    const DOCUMENTATION = "Documentation";
    await goToWMG(page);
    // verify logo
    expect(page.getByTestId("logo")).toBeVisible();

    // Help & Doc
    expect(await page.getByTestId("InputDropdown").textContent()).toBe(
      "Help & Documentation"
    );

    await page.getByTestId("InputDropdown").click();

    const popupPromise = page.waitForEvent("popup");
    await page.getByText(DOCUMENTATION, { exact: true }).click();
    const popup = await popupPromise;
    // Wait for new tab to load.
    await popup.waitForLoadState();
    expect(popup.url()).toContain(DOCUMENTATION);
    expect(
      popup.locator(
        getById("gene-expression--query-gene-expression-across-tissues")
      )
    ).toBeVisible();
  });
  test("Should take user to portal page on clicking on logo", async ({
    page,
  }) => {
    await tryUntil(
      async () => {
        await goToWMG(page);

        await page.getByTestId("logo").click();
        await page.waitForLoadState();

        await expect(
          page.getByTestId("laptop-with-cell-data-on-screen")
        ).toBeVisible();
      },
      { page }
    );
  });
});
