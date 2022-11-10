import { devices, expect, PlaywrightTestConfig } from "@playwright/test";
import { matchers } from "expect-playwright";
import fs from "fs";
import featureFlags from "./tests/common/featureFlags";

expect.extend(matchers);

const isHeadful =
  process.env.HEADFUL === "true" || process.env.HEADLESS === "false";

// 'github' for GitHub Actions CI to generate annotations, default otherwise
const PLAYWRIGHT_REPORTER = process.env.CI ? "github" : "list";

const VIEWPORT = {
  height: 1080,
  width: 1920,
};

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
  outputDir: "playwright-report/",

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
    acceptDownloads: true,
    /* Maximum time each action such as `click()` can take. Defaults to 0 (no limit). */
    actionTimeout: 0,
    headless: !isHeadful,
    ignoreHTTPSErrors: true,
    screenshot: "only-on-failure",
    storageState: { ...loginState(), ...featureFlags },
    /* Collect trace when retrying the failed test. See https://playwright.dev/docs/trace-viewer */
    trace: "retain-on-failure",
    video: {
      mode: "retain-on-failure",
      size: VIEWPORT,
    },
    viewport: VIEWPORT,
  },

  /* Opt out of parallel tests. */
  workers: 1,

  /* Run your local dev server before starting the tests */
  // webServer: {
  //   command: 'npm run start',
  //   port: 3000,
  // },
};

function loginState(): {
  cookies: Array<{
    name: string;

    value: string;

    /**
     * domain and path are required
     */
    domain: string;

    /**
     * domain and path are required
     */
    path: string;

    /**
     * Unix time in seconds.
     */
    expires: number;

    httpOnly: boolean;

    secure: boolean;

    /**
     * sameSite flag
     */
    sameSite: "Strict" | "Lax" | "None";
  }>;
} {
  let storageState = { cookies: [] };

  if (fs.existsSync("storageState.json")) {
    storageState = JSON.parse(fs.readFileSync("storageState.json", "utf-8"));
  }

  return storageState;
}

export default config;
