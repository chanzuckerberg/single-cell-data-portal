import { expect, test } from "@playwright/test";
import { goToPage } from "tests/utils/helpers";

const { describe } = test;

describe("Homepage", () => {
  test("renders the expected elements", async ({ page }) => {
    await goToPage(undefined, page);
    await expect(page.getByTestId("logo")).toHaveCount(2);
    await expect(page.getByTestId("collection-link")).toBeVisible();

    await Promise.all([
      page.waitForNavigation(),
      page.getByTestId("collection-link").click(),
    ]);

    expect(page.url()).toContain("collections");
  });
});
