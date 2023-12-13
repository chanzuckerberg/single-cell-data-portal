import { expect } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { goToPage } from "tests/utils/helpers";
import { TEST_URL } from "../../common/constants";
import { test } from "tests/common/test";

const { describe } = test;

describe("ToS and Privacy tests", () => {
  test("Should verify the ToS Page", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.TOS}`, page);

    await expect(page.getByText("Terms of Use").first()).toBeVisible();
    await expect(page.getByTestId("cellxgene-logo")).toBeVisible();
  });

  test("Should verify Privacy Page", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.PRIVACY}`, page);

    await expect(page.getByText("Privacy Policy").first()).toBeVisible();
    await expect(page.getByTestId("cellxgene-logo")).toBeVisible();
  });
});
