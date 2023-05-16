import { expect, Page, test } from "@playwright/test";
import { ROUTES } from "src/common/constants/routes";
import { goToPage, isDevStagingProd, tryUntil } from "tests/utils/helpers";
import { TEST_URL } from "../../common/constants";

const { describe, skip } = test;

describe("Cell Cards", () => {
  skip(!isDevStagingProd, "WMG BE API does not work locally or in rdev");

  describe("Landing Page", () => {
    test("All LandingPage components are present", async ({ page }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}`, page);
    });
    test("Cell type search bar filters properly and links to a CellCard", async ({
      page,
    }) => {
      await goToPage(`${TEST_URL}${ROUTES.CELL_CARDS}`, page);
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

async function waitForElement(page: Page, testId: string) {
  await tryUntil(
    async () => {
      await expect(page.getByTestId(testId)).not.toHaveCount(0);
    },
    { page }
  );
}

async function getElementAndClick(page: Page, testID: string) {
  await tryUntil(
    async () => {
      await page.getByTestId(testID).click();
    },
    { page }
  );
}
