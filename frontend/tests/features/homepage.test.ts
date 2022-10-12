import { expect, test } from "@playwright/test";
import { goToPage } from "tests/utils/helpers";
import { getTestID } from "../utils/selectors";

const { describe } = test;

describe("Homepage", () => {
  test("renders the expected elements", async ({ page }) => {
    await goToPage(undefined, page);
    await expect(page).toHaveSelectorCount(getTestID("logo"), 2);
    await expect(page).toHaveSelector(getTestID("collection-link"));

    await Promise.all([
      page.waitForNavigation(),
      page.click(getTestID("collection-link")),
    ]);

    expect(page.url()).toContain("collections");
  });
});
