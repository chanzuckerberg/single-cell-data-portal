import { test, chromium } from "@playwright/test";
import { SKIP_LOGIN } from "tests/common/constants";
import { COMMON_PLAYWRIGHT_CONTEXT } from "tests/common/context";
import featureFlags from "tests/common/featureFlags";
import { goToPage, login } from "tests/utils/helpers";

const { describe, skip } = test;

describe("global setup", () => {
  // conditional login based on SKIP_LOGIN
  test("login", async () => {
    skip(SKIP_LOGIN, "SKIP_LOGIN was set to true. Skipping log in.");

    const browser = await chromium.launch();

    const browserContext = await browser.newContext({
      ...COMMON_PLAYWRIGHT_CONTEXT,
      storageState: featureFlags,
    });

    const page = await browserContext.newPage();

    await goToPage(undefined, page);

    const apiUrl = await page.evaluate(() => {
      return document.body.getAttribute("data-api-url");
    });

    process.env.API_URL = apiUrl || "";

    try {
      // Skip login for tests that don't require it
      if (SKIP_LOGIN) {
        console.log("SKIP_LOGIN was set to true. Skipping log in.");
      } else {
        // One time auth
        console.log("Logging in...");

        await login(page);

        console.log(`Login success!`);

        await browserContext.close();
        await browser.close();
      }
    } catch (error) {
      console.error(error);
    }
  });

  test("set API_URL", async ({ page }) => {
    await goToPage(undefined, page);

    const apiUrl = await page.evaluate(() => {
      return document.body.getAttribute("data-api-url");
    });

    process.env.API_URL = apiUrl || "";
  });
});
