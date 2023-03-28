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
    await page.goto(`${TEST_URL}${ROUTES.WHERE_IS_MY_GENE}`);
    await page.waitForLoadState("networkidle");
  });

  test.only("Should verify main panel components", async ({ page }) => {
    await expect(page.locator(getTestID("add-tissue"))).toBeVisible();
    await expect(page.locator(getTestID("add-gene"))).toBeVisible();
    await expect(page.getByTestId("survey-alert-id")).toContainText(ALERT);
    await expect(page.getByText("quick survey")).toHaveAttribute(
      "href",
      SURVEY_LINK
    );
  });
});
