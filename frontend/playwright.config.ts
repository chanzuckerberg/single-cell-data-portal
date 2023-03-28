import {
  devices,
  expect,
  PlaywrightTestConfig,
  ReporterDescription,
} from "@playwright/test";
import { matchers } from "expect-playwright";
import fs from "fs";
import { LOGIN_STATE_FILENAME } from "tests/common/constants";
import { COMMON_PLAYWRIGHT_CONTEXT } from "tests/common/context";
import featureFlags from "./tests/common/featureFlags";

expect.extend(matchers);

// 'github' for GitHub Actions CI to generate annotations, default otherwise
const PLAYWRIGHT_REPORTER = process.env.CI
  ? ([["github"], ["line"], ["allure-playwright"]] as ReporterDescription[])
  : "list";

/**
 * See https://playwright.dev/docs/test-configuration.
 */
const config: PlaywrightTestConfig = {
  expect: {
    /**
     * Maximum time expect() should wait for the condition to be met.
     * For example in `await expect(locator).toHaveText();`
     */
    timeout: 10000,
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

  retries: process.env.CI ? 2 : 0,

  /* The base directory, relative to the config file, for snapshot files created with toMatchSnapshot and toHaveScreenshot. */
  snapshotDir: "./__snapshots__",

  testDir: "tests",

  /* Maximum time one test can run for. */
  timeout: 3 * 60 * 1000,

  use: {
    ...COMMON_PLAYWRIGHT_CONTEXT,
    storageState: getStorageState(),
  },

  /* Opt out of parallel tests. */
  workers: 1,

  /* Run your local dev server before starting the tests */
  // webServer: {
  //   command: 'npm run start',
  //   port: 3000,
  // },
};

function getStorageState(): {
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
  origins: Array<{
    origin: string;

    localStorage: Array<{
      name: string;

      value: string;
    }>;
  }>;
} {
  const storageState = featureFlags;

  if (fs.existsSync(LOGIN_STATE_FILENAME)) {
    const loginState = JSON.parse(
      fs.readFileSync(LOGIN_STATE_FILENAME, "utf-8")
    );

    // Merge loginState with featureFlags
    storageState.cookies = storageState.cookies.concat(loginState.cookies);
    storageState.origins = storageState.origins.concat(loginState.origins);
  }

  return storageState;
}

export default config;
