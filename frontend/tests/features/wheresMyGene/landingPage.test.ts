import { expect, Page } from "@playwright/test";
import { TEST_URL } from "../../common/constants";
import { ROUTES } from "src/common/constants/routes";
import { ADD_GENE_BTN } from "tests/common/constants";
import { getById } from "tests/utils/selectors";
import { goToPage, tryUntil } from "tests/utils/helpers";
import { test } from "tests/common/test";
import {
  waitForHeatmapToRender,
  WMG_WITH_SEEDED_GENES,
} from "tests/utils/wmgUtils";

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

describe("Tests for Citation, Share and Notifications", () => {
  test("Citation Button - check visibility, notification and disabled state", async ({
    page,
  }) => {
    /**
     * Check Citation button is enabled on landing page without genes
     */
    await goToWMG(page);
    await expect(page.getByTestId(LEGEND_WRAPPER)).toContainText("Citation");
    const citationButton = page.getByTestId("citation-button");
    await expect(citationButton).toBeEnabled();
    /**
     * Check Citation button is disabled after click and notification is shown
     */
    await citationButton.click();
    await expect(citationButton).toBeDisabled();
    const notification = page.getByTestId("notification");
    await expect(notification).toBeVisible();
    await expect(notification).toContainText("Citation copied to clipboard");
    /**
     * Check that notification is removed after timeout and button is enabled
     */
    await tryUntil(
      async () => {
        await expect(notification).toBeHidden();
        await expect(citationButton).toBeEnabled();
      },
      { page }
    );
  });
  test("Share Button - check visibility, notification and disabled state", async ({
    page,
  }) => {
    /**
     * Check Share button is disabled on landing page without genes
     */
    await goToWMG(page);
    await expect(page.getByTestId(LEGEND_WRAPPER)).toContainText("Share");
    const shareButton = page.getByTestId("share-button");
    await expect(shareButton).toBeDisabled();
    /**
     * Check Share button is enabled on landing page with genes
     */
    await goToPage(WMG_WITH_SEEDED_GENES.URL, page);
    await waitForHeatmapToRender(page);
    await expect(shareButton).toBeEnabled();
    /**
     * Check Share button is disabled after clicking and notification is shown
     */
    await shareButton.click();
    await expect(shareButton).toBeDisabled();
    const notification = page.getByTestId("notification");
    await expect(notification).toBeVisible();
    await expect(notification).toContainText("Share URL copied to clipboard");
    /**
     * Check that notification is removed after timeout and button is enabled
     */
    await tryUntil(
      async () => {
        await expect(notification).toBeHidden();
        await expect(shareButton).toBeEnabled();
      },
      { page }
    );
    /**
     * Check that notifications stack and are removed in order for both buttons
     */
    const citationButton = page.getByTestId("citation-button");
    await citationButton.click();
    await shareButton.click();
    await expect(notification).toHaveCount(2);
    await tryUntil(
      async () => {
        await expect(notification).toHaveCount(0);
      },
      { page }
    );
  });
});
