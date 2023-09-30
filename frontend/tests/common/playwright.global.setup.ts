import { test, chromium } from "@playwright/test";
import { SKIP_LOGIN } from "tests/common/constants";
import { COMMON_PLAYWRIGHT_CONTEXT } from "tests/common/context";
import { getFeatureFlags } from "tests/common/featureFlags";
import { goToPage, isDevStagingRdev, login } from "tests/utils/helpers";

const { describe, skip } = test;

describe("global setup", () => {
  test("login @loggedIn", async () => {
    skip(SKIP_LOGIN, "SKIP_LOGIN was set to true. Skipping log in.");
    skip(
      !isDevStagingRdev,
      `Skipping login for environment that is not "dev", "rdev', or "staging"`
    );

    const browser = await chromium.launch();
    const browserContext = await browser.newContext({
      ...COMMON_PLAYWRIGHT_CONTEXT,
      storageState: getFeatureFlags({ curator: true }),
    });

    const page = await browserContext.newPage();

    await goToPage(undefined, page);

    // One time auth
    await login(page);

    await browserContext.close();
    await browser.close();
  });
});
