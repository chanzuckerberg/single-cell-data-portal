import { expect, Page, test } from "@playwright/test";
import { TEST_URL } from "../../common/constants";
import { ROUTES } from "src/common/constants/routes";

const { describe } = test;
const ALERT =
  "We would appreciate your feedback, please fill out a quick survey";

const SURVEY_LINK = "https://airtable.com/shrLwepDSEX1HI6bo";
const ADD_TISSUE = "add-tissue";
const ADD_GENES = "add-genes";
const EXPLORE_GENE_EXPRESSION = "explore-gene-expression";
const LEGEND_WRAPPER = "legend-wrapper";
const DOT_SIZES: Record<number, number> = {
  0: 4,
  1: 9,
  2: 12,
  3: 14,
  4: 16,
};
async function checkDotSize(index: number, page: Page) {
  expect(
    page
      .locator('[id="expressed-in-cells-dots"]')
      .locator("span")
      .nth(0)
      .getAttribute("size")
  ).toBe(DOT_SIZES[index]);
}
describe("Tests for Gene Expression page", () => {
  test.beforeEach(async ({ page }) => {
    await Promise.all([
      page.waitForResponse(
        (resp: { url: () => string | string[]; status: () => number }) =>
          resp.url().includes("/wmg/v1/filters") && resp.status() === 200
      ),
      page.goto(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`),
    ]);
  });

  test("Should verify main panel components", async ({ page }) => {
    // +Tissue button
    await expect(page.getByTestId(ADD_TISSUE).first()).toBeVisible();
    // +Gene button
    await expect(page.getByTestId(ADD_GENES).first()).toBeVisible();
    // survey alert
    await expect(page.getByTestId("survey-alert-id")).toContainText(ALERT);
    await expect(page.getByText("quick survey")).toHaveAttribute(
      "href",
      SURVEY_LINK
    );
    // default organism filter
    await expect(
      page.getByTestId("add-organism").locator("span")
    ).toContainText("Homo sapiens");

    // STEP 1 column
    await expect(page.getByTestId("column-one")).toBeVisible();
    await expect(page.getByTestId(ADD_TISSUE).nth(1)).toContainText("STEP 1");
    await expect(page.getByTestId(ADD_TISSUE).nth(1)).toContainText(
      "Add Tissues"
    );

    // STEP 2 column
    await expect(page.getByTestId("column-two")).toBeVisible();
    await expect(page.getByTestId(ADD_GENES).nth(1)).toContainText("STEP 2");
    await expect(page.getByTestId(ADD_GENES).nth(1)).toContainText("Add Genes");

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
      page.locator('[data-nimg="intrinsic"]').locator("svg")
    ).toBeVisible();

    // Gene expression in cells
    await expect(page.getByTestId(LEGEND_WRAPPER)).toContainText(
      "Expressed in Cells (%)"
    );
    await expect(page.locator('[id="expressed-in-cells-dots"]')).toBeVisible();

    Object.keys(DOT_SIZES).forEach(async (index) => {
      await checkDotSize(Number(index), page);
    });
  });
});
