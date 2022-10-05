import type { PlaywrightTestConfig } from "@playwright/test";
import { devices, expect } from "@playwright/test";
import { matchers } from "expect-playwright";
import featureFlags from "./tests/common/featureFlags";

expect.extend(matchers);

const isHeadful =
  process.env.HEADFUL === "true" || process.env.HEADLESS === "false";

const DEFAULT_CONTEXT_CONFIG = {
  acceptDownloads: true,
  args: ["--ignore-certificate-errors", "--ignore-ssl-errors"],
  headless: !isHeadful,
  ignoreHTTPSErrors: true,
  storageState: featureFlags,
};

// 'github' for GitHub Actions CI to generate annotations, default otherwise
const PLAYWRIGHT_REPORTER = process.env.CI ? "github" : "list";

/**
 * See https://playwright.dev/docs/test-configuration.
 */
const config: PlaywrightTestConfig = {
  expect: {
    /**
     * Maximum time expect() should wait for the condition to be met.
     * For example in `await expect(locator).toHaveText();`
     */
    timeout: 5000,
  },

  /* Fail the build on CI if you accidentally left test.only in the source code. */
  forbidOnly: !!process.env.CI,

  /* Run tests in files in parallel */
  fullyParallel: true,

  globalSetup: require.resolve("./playwright-globalSetup"),

  /* Folder for test artifacts such as screenshots, videos, traces, etc. */
  outputDir: "test-results/",

  /* Configure projects for major browsers */
  projects: [
    {
      name: "chromium",
      use: {
        ...devices["Desktop Chrome"],
        /**
         * (thuang): Add `czi-checker`, so Plausible will ignore it.
         * NOTE: This changes all browsers to use Desktop Chrome UA, so please look
         * out for bugs that could be caused by this.
         * https://github.com/matomo-org/device-detector/blob/master/regexes/bots.yml#L2762
         */
        userAgent: devices["Desktop Chrome"].userAgent + " czi-checker",
      },
    },
  ],

  /* Reporter to use. See https://playwright.dev/docs/test-reporters */
  reporter: PLAYWRIGHT_REPORTER,

  retries: 2,

  /* The base directory, relative to the config file, for snapshot files created with toMatchSnapshot and toHaveScreenshot. */
  snapshotDir: "./__snapshots__",

  testDir: "tests",

  /* Maximum time one test can run for. */
  timeout: 3 * 60 * 1000,

  use: {
    ...DEFAULT_CONTEXT_CONFIG,
    /* Maximum time each action such as `click()` can take. Defaults to 0 (no limit). */
    actionTimeout: 0,
    /* Base URL to use in actions like `await page.goto('/')`. */
    // baseURL: 'http://localhost:3000',
    screenshot: "only-on-failure",
    /* Collect trace when retrying the failed test. See https://playwright.dev/docs/trace-viewer */
    trace: "on-first-retry",
    video: { mode: "retain-on-failure", size: { height: 1080, width: 1920 } },
  },

  /* Opt out of parallel tests on CI. */
  workers: process.env.CI ? 1 : undefined,

  /* Run your local dev server before starting the tests */
  // webServer: {
  //   command: 'npm run start',
  //   port: 3000,
  // },
};

export default config;
