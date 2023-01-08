import { chromium } from "@playwright/test";
import { SKIP_LOGIN, TEST_URL } from "tests/common/constants";
import { COMMON_PLAYWRIGHT_CONTEXT } from "tests/common/context";
import featureFlags from "tests/common/featureFlags";
import { login } from "tests/utils/helpers";

module.exports = async () => {
  try {
    // Skip login for tests that don't require it
    if (SKIP_LOGIN) {
      console.log("SKIP_LOGIN was set to true. Skipping log in.");
    } else {
      // One time auth
      const browser = await chromium.launch();
      const browserContext = await browser.newContext({
        ...COMMON_PLAYWRIGHT_CONTEXT,
        storageState: featureFlags,
      });
      const page = await browserContext.newPage();
      await page.goto(TEST_URL, { timeout: 60000 });

      console.log("Logging in...");

      await login(page);

      console.log(`Login success!`);

      await browserContext.close();
      await browser.close();
    }
  } catch (error) {
    console.error(error);
  }
};
