import { test, chromium } from "@playwright/test";
import { SKIP_LOGIN } from "tests/common/constants";
import { COMMON_PLAYWRIGHT_CONTEXT } from "tests/common/context";
import { getFeatureFlags } from "tests/common/featureFlags";
import { goToPage, isDevStaging, login } from "tests/utils/helpers";

const { describe, skip } = test;

describe("global setup", () => {
  test("login", async () => {
    skip(SKIP_LOGIN, "SKIP_LOGIN was set to true. Skipping log in.");
    skip(
      !isDevStaging,
      `Skipping login for environment that is not "dev" or "staging"`
    );

    const browser = await chromium.launch();
    const browserContext = await browser.newContext({
      ...COMMON_PLAYWRIGHT_CONTEXT,
      storageState: getFeatureFlags({ curator: true }),
    });

    const page = await browserContext.newPage();

    await goToPage(undefined, page);

    const apiUrl = await page.evaluate(() => {
      return document.body.getAttribute("data-api-url");
    });

    process.env.API_URL = apiUrl || "";

    // One time auth
    await login(page);

    await browserContext.close();
    await browser.close();
  });

  test("set API_URL", async ({ page }) => {
    await goToPage(undefined, page);

    const apiUrl = await page.evaluate(() => {
      return document.body.getAttribute("data-api-url");
    });

    process.env.API_URL = apiUrl || "";
  });
});
