import { expect, test } from "@playwright/test";
import { TEST_URL } from "../../common/constants";
import { ROUTES } from "src/common/constants/routes";
import { getTestID } from "tests/utils/selectors";

const { describe } = test;
const ALERT =
  "We would appreciate your feedback, please fill out a quick survey";

const SURVEY_LINK = "https://airtable.com/shrLwepDSEX1HI6bo";

describe("Tests for Gene Expression page", () => {
  test.beforeEach(async ({ page }) => {
    // await page.goto(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`);
    // await page.waitForLoadState("networkidle");
    await Promise.all([
      page.waitForResponse(
        (resp: { url: () => string | string[]; status: () => number }) =>
          resp.url().includes("/wmg/v1/filters") && resp.status() === 200
      ),
      page.goto(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`),
    ]);
  });

  test.only("Should verify main panel components", async ({ page }) => {
    // +Tissue button
    await expect(page.getByTestId("add-tissue")).toBeVisible();
    // +Gene button
    await expect(page.getByTestId("add-gene")).toBeVisible();
    // survey alert
    await expect(page.getByTestId("survey-alert-id")).toContainText(ALERT);
    await expect(page.getByText("quick survey")).toHaveAttribute(
      "href",
      SURVEY_LINK
    );
    // default organism filter
    await expect(page.getByTestId("add-organism").locator("span")).toContainText(
      "Homo sapiens"
    );

    // add tissues column
    await expect(page.getByTestId("column-one")).toBeVisible();
    await expect(page.getByTestId("add-tissue")).toContainText("Step 1");
    await expect(page.getByTestId("add-tissue")).toContainText("Add Tissues");
  });
  test("Should verify organism filter", async ({ page }) => {
    // +Tissue button
    await expect(page.getByTestId("add-tissue")).toBeVisible();
    // +Gene button
    await expect(page.locator(getTestID("add-gene"))).toBeVisible();
    // survey alert
    await expect(page.getByTestId("survey-alert-id")).toContainText(ALERT);
    await expect(page.getByText("quick survey")).toHaveAttribute(
      "href",
      SURVEY_LINK
    );
    // default organism filter
    expect(page.locator('[data-test-id="addorganism"]')).toContainText(
      "Homo sapiens"
    );
  });
});
