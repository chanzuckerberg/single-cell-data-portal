import { expect, test } from "@playwright/test";
import { TEST_URL } from "tests/common/constants";
import { goToPage, tryUntil } from "tests/utils/helpers";
import { getTestID } from "tests/utils/selectors";

const { describe } = test;

describe("component Chart", () => {
  test("the expected elements are rendered", async ({ page }) => {
    await goToPage(`${TEST_URL}/tests/components/Chart`, page);

    await tryUntil(
      async () => {
        await expect(page.getByTestId("chart")).toBeVisible();
      },
      { page }
    );

    /**
     * (thuang): Since macos only generates screenshot files with `darwin` suffix,
     * which don't work for CI/CD Linux environment, we need to manually duplicate
     * the generated darwin files and rename them to linux.
     * See: https://github.com/microsoft/playwright/discussions/12729
     */
    await expect(page.locator(getTestID("chart"))).toHaveScreenshot();
    await expect(page.locator(getTestID("canvas-chart"))).toHaveScreenshot();
  });
});
