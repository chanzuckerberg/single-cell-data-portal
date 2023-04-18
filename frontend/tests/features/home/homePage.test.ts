import { Page, expect, test } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "tests/common/constants";
import { goToPage, tryUntil } from "tests/utils/helpers";

const { describe } = test;
const COLLECTION_LINK = "collection-link";
describe("Homepage", () => {
  test("Should verify the expected elements are rendered", async ({ page }) => {
    await goToPage(undefined, page);
    expect(await page.getByTestId("logo").count()).toBe(2);
    await expect(page.getByTestId(COLLECTION_LINK)).toBeVisible();

    await isPageScrollableToSeeSiteMap(page);

    // check collection link works
    await page.getByTestId(COLLECTION_LINK).click();
    await page.waitForURL("**" + ROUTES.COLLECTIONS);

    await isGlobalLayoutWrapperScrollable(page);
  });
  test("Should render the ToS Page", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.TOS}`, page);

    await tryUntil(
      async () => {
        await expect(
          page.locator("h1").getByText("Terms of Use")
        ).toBeVisible();
      },
      { page }
    );

    await expect(page.getByTestId("cellxgene-logo")).toBeVisible();
    await isGlobalLayoutWrapperScrollable(page);
  });

  test("Should render the Privacy Page", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.PRIVACY}`, page);

    await tryUntil(
      async () => {
        await expect(
          page.locator("h1").getByText("Privacy Policy")
        ).toBeVisible();
      },
      { page }
    );

    await expect(page.getByTestId("cellxgene-logo")).toBeVisible();

    await isGlobalLayoutWrapperScrollable(page);
  });

  test("Should render Sitemap", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.SITEMAP}`, page);

    await tryUntil(
      async () => {
        await expect(page.locator("h1").getByText("Sitemap")).toBeVisible();
      },
      { page }
    );

    await isPageScrollableToSeeSiteMap(page);
  });
  test("Should render collections page", async ({ page }) => {
    await goToPage(`${TEST_URL}${ROUTES.COLLECTIONS}`, page);

    await isGlobalLayoutWrapperScrollable(page);
  });

  test("Should render documentation page", async ({ page, context }) => {
    await goToPage(undefined, page);

    await page.getByText("Help & Documentation").click();

    // New tab
    const newPage = await context.waitForEvent("page");

    await tryUntil(
      async () => {
        await expect(
          newPage.getByText("Welcome to CZ CELLxGENE")
        ).toBeVisible();
      },
      { page: newPage }
    );

    await isGlobalLayoutWrapperScrollable(newPage);
  });
});

async function isPageScrollableToSeeSiteMap(page: Page) {
  const SitemapLink = page.locator("a", { hasText: "Sitemap" });

  await expect(SitemapLink).not.toBeInViewport();
  await page.mouse.wheel(0, 20000);
  await expect(SitemapLink).toBeInViewport();
}

/**
 * (thuang): Some pages don't have sitemap in the footer, so this is the alternative
 * to check if the page is scrollable
 */
async function isGlobalLayoutWrapperScrollable(page: Page) {
  const wrapper = page.getByTestId("global-layout-wrapper");

  expect(
    await wrapper.evaluate((e) => {
      return e.scrollTop === 0;
    })
  ).toBeTruthy();

  await tryUntil(
    async () => {
      // (thuang): Click is needed to focus the scrollable element
      await wrapper.click();
      await page.mouse.wheel(0, 20000);

      expect(
        await wrapper.evaluate((e) => {
          return e.scrollTop > e.clientHeight;
        })
      ).toBeTruthy();
    },
    { page }
  );
}
