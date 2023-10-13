import { expect, Page } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { TEST_URL } from "tests/common/constants";
import { goToPage, tryUntil } from "tests/utils/helpers";
import {
  LANDING_PAGE_FALLBACK_CELLS_HERO_NUM,
  LANDING_PAGE_FALLBACK_CELLTYPES_HERO_NUM,
  LANDING_PAGE_FALLBACK_DATASETS_HERO_NUM,
  LANDING_PAGE_CELLS_HERO_NUM_ID,
  LANDING_PAGE_CELLTYPES_HERO_NUM_ID,
  LANDING_PAGE_DATASETS_HERO_NUM_ID,
} from "src/views/Landing/constants";
import { test } from "tests/common/test";

const { describe } = test;
const COLLECTIONS_LINK_ID = "collections-link";
const SCROLL_Y_PX = 999999;

describe("Homepage", () => {
  test("Should verify the expected elements are rendered", async ({ page }) => {
    await goToPage(undefined, page);
    await tryUntil(
      async () => {
        expect(await page.getByTestId("logo").count()).toBe(2);
      },
      { page }
    );
    await expect(page.getByTestId(COLLECTIONS_LINK_ID)).toBeVisible();

    await isPageScrollableToSeeSiteMap(page);

    // check collection link works
    await page.getByTestId(COLLECTIONS_LINK_ID).click();
    await page.waitForURL("**" + ROUTES.COLLECTIONS);

    await isGlobalLayoutWrapperScrollable(page);
  });
  test("Hero numbers are rendered", async ({ page }) => {
    await goToPage(undefined, page);

    await tryUntil(
      async () => {
        const cellsHeroNum = await page
          .getByTestId(LANDING_PAGE_CELLS_HERO_NUM_ID)
          .innerText();
        const cellTypesHeroNum = await page
          .getByTestId(LANDING_PAGE_CELLTYPES_HERO_NUM_ID)
          .innerText();
        const datasetsHeroNum = await page
          .getByTestId(LANDING_PAGE_DATASETS_HERO_NUM_ID)
          .innerText();

        expect(cellsHeroNum).not.toEqual(LANDING_PAGE_FALLBACK_CELLS_HERO_NUM);
        expect(cellTypesHeroNum).not.toEqual(
          LANDING_PAGE_FALLBACK_CELLTYPES_HERO_NUM
        );
        expect(datasetsHeroNum).not.toEqual(
          LANDING_PAGE_FALLBACK_DATASETS_HERO_NUM
        );
      },
      { page }
    );
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
  const Html = page.locator("html");
  const Body = page.locator("body");

  const bodyHeight = await Body.evaluate((e) => {
    return e.clientHeight;
  });
  const viewportHeight = await Html.evaluate((e) => {
    return e.clientHeight;
  });

  const scrollTopTarget = bodyHeight - viewportHeight - 50;

  await expect(SitemapLink).not.toBeInViewport();

  /**
   * (thuang): Scroll to the bottom iteratively, because Firefox somehow has a
   * max scroll height limit, so we can't just scroll to the bottom in one go
   */
  await tryUntil(
    async () => {
      await page.mouse.wheel(0, SCROLL_Y_PX);

      const scrollTop = await Html.evaluate(async (e) => {
        return e.scrollTop;
      });

      expect(scrollTop > scrollTopTarget).toBeTruthy();
    },
    { page }
  );

  await expect(SitemapLink).toBeInViewport();
}

/**
 * (thuang): Some pages don't have sitemap in the footer, so this is the alternative
 * to check if the page is scrollable
 */
async function isGlobalLayoutWrapperScrollable(page: Page) {
  if (
    await page
      .locator("main")
      .evaluate((e) => e.scrollHeight <= window.innerHeight)
  ) {
    /**
     * (thuang): Set the viewport size to a smaller size, so that we can scroll
     */
    await page.setViewportSize({ height: 300, width: 600 });
  }

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
      await page.mouse.wheel(0, SCROLL_Y_PX);

      expect(
        await wrapper.evaluate((e) => {
          return e.scrollTop > 0;
        })
      ).toBeTruthy();
    },
    { page }
  );
}
