import { test, Page, expect, Locator } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { goToPage, tryUntil } from "tests/utils/helpers";
import { TEST_URL } from "../../common/constants";
import { LANDING_PAGE_HEADER } from "src/views/CellCards/components/LandingPage";
import {
  CELL_CARD_SEARCH_BAR,
  CELL_CARD_SEARCH_BAR_TEXT_INPUT,
} from "src/views/CellCards/components/CellCardSearchBar";
import {
  CELL_CARD_CL_DESCRIPTION,
  CELL_CARD_GPT_DESCRIPTION,
  CELL_CARD_GPT_TOOLTIP_LINK,
} from "src/views/CellCards/components/CellCard/components/Description";
import {
  CELL_CARD_HEADER_NAME,
  CELL_CARD_HEADER_TAG,
} from "src/views/CellCards/components/CellCard";

const { describe } = test;

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
    test("All cell card components are present", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page);
      await isElementVisible(page, CELL_CARD_HEADER_NAME);
      await isElementVisible(page, CELL_CARD_HEADER_TAG);
      await isElementVisible(page, CELL_CARD_CL_DESCRIPTION);
      await isElementVisible(page, CELL_CARD_GPT_DESCRIPTION);
      await isElementVisible(page, CELL_CARD_GPT_TOOLTIP_LINK);
      await isElementVisible(page, CELL_CARD_SEARCH_BAR);
      const headerName = page.getByTestId(CELL_CARD_HEADER_NAME);
      const headerNameText = await headerName.textContent();
      expect(headerNameText).toBe("Neuron");
    });
    test("Cell card GPT description tooltip displays disclaimer", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page);
      await isElementVisible(page, CELL_CARD_GPT_TOOLTIP_LINK);
      await page.getByTestId(CELL_CARD_GPT_TOOLTIP_LINK).hover();
      // check role tooltip is visible
      const tooltipLocator = page.getByRole("tooltip");
      await tooltipLocator.waitFor({ timeout: 5000 });
      const tooltipLocatorVisible = await tooltipLocator.isVisible();
      expect(tooltipLocatorVisible).toBe(true);
      // check that tooltip contains disclaimer
      const tooltipText = await tooltipLocator.textContent();
      expect(tooltipText).toContain(
        `This summary on "neuron" was generated with ChatGPT, powered by the GPT3.5 Turbo model. Keep in mind that ChatGPT may occasionally present information that is not entirely accurate. For transparency, the prompts used to generate this summary are shared below.`
      );
    });
    test("Cell type search bar filters properly and links to a CellCard", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000540`, page);
      const element = page.getByTestId(CELL_CARD_SEARCH_BAR_TEXT_INPUT);
      await element.waitFor();
      await element.click();
      // get number of elements with role option in dropdown
      const numOptionsBefore = await countLocator(page.getByRole("option"));
      // type in search bar
      await element.type("acinar cell");
      // get number of elements with role option in dropdown
      const numOptionsAfter = await countLocator(page.getByRole("option"));
      // check that number of elements with role option in dropdown has decreased
      expect(numOptionsAfter).toBeLessThan(numOptionsBefore);
      // check that the first element in the dropdown is the one we searched for
      const firstOption = (await page.getByRole("option").elementHandles())[0];
      const firstOptionText = await firstOption?.textContent();
      expect(firstOptionText).toBe("acinar cell");
      // click on first element in dropdown
      await firstOption?.click();
      // check that the url has changed to the correct cell card
      await page.waitForURL(`${TEST_URL}${ROUTES.CELL_CARDS}/CL_0000622`);
    });
  });
});

async function isElementVisible(page: Page, testId: string) {
  await tryUntil(
    async () => {
      const element = page.getByTestId(testId);
      await element.waitFor({ timeout: 5000 });
      const isVisible = await element.isVisible();
      expect(isVisible).toBe(true);
    },
    { page }
  );
}

async function countLocator(locator: Locator) {
  return (await locator.elementHandles()).length;
}
