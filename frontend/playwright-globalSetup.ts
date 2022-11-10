import { chromium } from "@playwright/test";
import { TEST_URL } from "tests/common/constants";
import featureFlags from "tests/common/featureFlags";
import { login } from "tests/utils/helpers";

module.exports = async () => {
  try {
    // One time auth
    const browser = await chromium.launch();
    const browserContext = await browser.newContext({
      storageState: featureFlags,
    });
    const page = await browserContext.newPage();
    await page.goto(TEST_URL);

    console.log("Logging in...");

    await login(page);

    console.log(`Login success!`);

    await browserContext.close();
    await browser.close();
  } catch (error) {
    console.error(error);
  }
};
