import { expect, test } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { goToPage } from "tests/utils/helpers";
import { TEST_URL } from "../common/constants";
import { getTestID, getText } from "../utils/selectors";

const { describe } = test;

describe("ToS and Privacy", () => {
  test("renders the the ToS Page", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.TOS}`, page);

    await expect(page).toHaveSelector(getText("Terms of Use"));
    await expect(page).toHaveSelector(getTestID("cellxgene-logo"));
  });

  test("renders the the Privacy Page", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.PRIVACY}`, page);

    await expect(page).toHaveSelector(getText("Privacy Policy"));
    await expect(page).toHaveSelector(getTestID("cellxgene-logo"));
  });
});
