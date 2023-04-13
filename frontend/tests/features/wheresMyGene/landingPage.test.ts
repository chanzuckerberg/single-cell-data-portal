import { expect, Page, test } from "@playwright/test";
import { TEST_URL } from "../../common/constants";
import { ROUTES } from "src/common/constants/routes";
import {
  ADD_GENE_BTN,
  ADD_GENE_LBL,
  ADD_TISSUE_LBL,
} from "tests/utils/constants";

const ALERT =
  "We would appreciate your feedback, please fill out a quick survey";

const SURVEY_LINK = "https://airtable.com/shrLwepDSEX1HI6bo";
const EXPLORE_GENE_EXPRESSION = "explore-gene-expression";
const LEGEND_WRAPPER = "legend-wrapper";
const DOT_SIZES = ["4", "9", "12", "14", "16"];

function goToWMG(page: Page) {
  return Promise.all([
    page.waitForResponse(
      (resp: { url: () => string | string[]; status: () => number }) =>
        resp.url().includes("/wmg/v1/filters") && resp.status() === 200
    ),
    page.goto(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`),
  ]);
}
test.describe("Tests for Gene Expression page", () => {
  test("Should verify main panel components", async ({ page }) => {
    await goToWMG(page);
    // +Tissue button
    await expect(page.getByTestId(ADD_GENE_BTN)).toBeVisible();
    // +Gene button
    await expect(page.getByTestId(ADD_GENE_BTN)).toBeVisible();
    // survey alert
    await expect(page.getByTestId("survey-alert-id")).toContainText(ALERT);
    await expect(page.getByText("quick survey")).toHaveAttribute(
      "href",
      SURVEY_LINK
    );
    // default organism filter
    await expect(page.getByTestId("add-organism")).toContainText(
      "Homo sapiens"
    );

    // STEP 1 column
    await expect(page.getByTestId("column-one")).toBeVisible();
    await expect(page.getByTestId(ADD_TISSUE_LBL)).toContainText("STEP 1");
    await expect(page.getByTestId(ADD_TISSUE_LBL)).toContainText("Add Tissues");

    // STEP 2 column
    await expect(page.getByTestId("column-two")).toBeVisible();
    await expect(page.getByTestId(ADD_GENE_LBL)).toContainText("STEP 2");
    await expect(page.getByTestId(ADD_GENE_LBL)).toContainText("Add Genes");

    // STEP 3 column
    await expect(page.getByTestId(EXPLORE_GENE_EXPRESSION)).toContainText(
      "STEP 3"
    );
    await expect(page.getByTestId(EXPLORE_GENE_EXPRESSION)).toContainText(
      "Explore Gene Expression"
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
});
