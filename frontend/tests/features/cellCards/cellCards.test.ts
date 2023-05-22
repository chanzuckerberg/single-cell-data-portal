import { test, Page, expect, Locator } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { goToPage, tryUntil } from "tests/utils/helpers";
import { TEST_URL } from "../../common/constants";

const { describe } = test;

const LANDING_PAGE_HEADER = "landing-page-header";
const CELL_CARD_SEARCH_BAR = "cell-card-search-bar";
const CELL_CARD_SEARCH_BAR_TEXT_INPUT = "cell-card-search-bar-text-input";

describe("Cell Cards", () => {
  describe("Landing Page", () => {
    test("All LandingPage components are present", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}`, page);
      await isElementVisible(page, LANDING_PAGE_HEADER);
      await isElementVisible(page, CELL_CARD_SEARCH_BAR);
    });
    test("Cell type search bar filters properly and links to a CellCard", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}`, page);
      const element = page.getByTestId(CELL_CARD_SEARCH_BAR_TEXT_INPUT);
      await element.waitFor();
      await element.click();
      // get number of elements with role option in dropdown
      const numOptionsBefore = await countLocator(page.getByRole("option"));
      // type in search bar
      await element.type("neuron");
      // get number of elements with role option in dropdown
      const numOptionsAfter = await countLocator(page.getByRole("option"));
      // check that number of elements with role option in dropdown has decreased
      expect(numOptionsAfter).toBeLessThan(numOptionsBefore);
      // check that the first element in the dropdown is the one we searched for
      const firstOption = (await page.getByRole("option").elementHandles())[0];
      const firstOptionText = await firstOption?.textContent();
      expect(firstOptionText).toBe("neuron");
      // click on first element in dropdown
      await firstOption?.click();
      // check that the url has changed to the correct cell card
      await page.waitForURL(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`);
    });
    test("Cell type search bar keyboard input works properly", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}`, page);
      const element = page.getByTestId(CELL_CARD_SEARCH_BAR_TEXT_INPUT);
      await element.waitFor();
      await element.click();
      await element.type("neuron");
      // input down arrow key
      await element.press("ArrowDown");
      // input enter
      await element.press("Enter");
      // check that the url has changed to the correct cell card after browser finishes navigating
      await page.waitForURL(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`);
    });
  });

  describe("Cell Card", () => {
    test("All CellCard components are present", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}`, page);
    });
    test("Cell type search bar filters properly and links to a CellCard", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}`, page);
    });
  });
});

async function isElementVisible(page: Page, testId: string) {
  await tryUntil(
    async () => {
      const element = page.getByTestId(testId);
      const isVisible = await element.isVisible();
      expect(isVisible).toBe(true);
    },
    { page }
  );
}

// async function getElementAndClick(page: Page, testID: string) {
//   await tryUntil(
//     async () => {
//       await page.getByTestId(testID).click();
//     },
//     { page }
//   );
// }

// (alec) use this instead of locator.count() to make sure that the element is actually present
// when counting
async function countLocator(locator: Locator) {
  return (await locator.elementHandles()).length;
}
