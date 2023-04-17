import { expect, test } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "tests/common/constants";
import { goToPage, scrollToPageBottom } from "tests/utils/helpers";

const { describe } = test;
const COLLECTION_LINK = "collection-link";
describe("Homepage", () => {
  test("Should verify the expected elements are rendered", async ({ page }) => {
    await goToPage(undefined, page);
    expect(await page.getByTestId("logo").count()).toBe(2);
    await expect(page.getByTestId(COLLECTION_LINK)).toBeVisible();

    // check page scroll
    expect(await scrollToPageBottom(page)).toBeTruthy;

    // check collection link works
    await Promise.all([
      page.waitForNavigation(),
      page.getByTestId(COLLECTION_LINK).click(),
    ]);
    expect(page.url()).toContain("collection");
    // check page scroll
    expect(await scrollToPageBottom(page)).toBeTruthy;
  });
  test("Should render the ToS Page", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.TOS}`, page);

    await expect(page.locator("h1").getByText("Terms of Use")).toBeVisible();
    await expect(page.getByTestId("cellxgene-logo")).toBeVisible();
    // check page scroll
    expect(await scrollToPageBottom(page)).toBeTruthy;
  });

  test("Should render the Privacy Page", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.PRIVACY}`, page);

    await expect(page.locator("h1").getByText("Privacy Policy")).toBeVisible();
    await expect(page.getByTestId("cellxgene-logo")).toBeVisible();
    // check page scroll
    expect(await scrollToPageBottom(page)).toBeTruthy;
  });

  test("Should render Sitemap", async ({ page }) => {
    await goToPage(undefined, page);

    await page.getByText("Sitemap").click();
    await expect(page.locator("h1").getByText("Sitemap")).toBeVisible();
    // check page scroll
    expect(await scrollToPageBottom(page)).toBeTruthy;
  });
  test("Should render collections page", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.COLLECTIONS}`, page);
    // check page scroll
    expect(await scrollToPageBottom(page)).toBeTruthy;
  });

  test("Should render documentation page", async ({ page, context }) => {
    await goToPage(undefined, page);

    const [newPage] = await Promise.all([
      context.waitForEvent("page"),
      page.getByText("Help & Documentation").click(), // Opens a new tab
    ]);
    await expect(
      newPage.locator('[id="welcome-to-cz-cellxgene"]')
    ).toBeVisible();
    // check page scroll
    expect(await scrollToPageBottom(page)).toBeTruthy;
  });
  //todo: tests for gx to be added to the gx main panel tests
});
