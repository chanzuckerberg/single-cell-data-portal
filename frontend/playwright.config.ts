import {
  devices,
  expect,
  PlaywrightTestConfig,
  ReporterDescription,
  defineConfig,
} from "@playwright/test";
import { matchers } from "expect-playwright";
import fs from "fs";
import { LOGIN_STATE_FILENAME } from "tests/common/constants";
import { COMMON_PLAYWRIGHT_CONTEXT } from "tests/common/context";
import { getFeatureFlags } from "tests/common/featureFlags";

const CICD_MAX_FAILURE = 2;

expect.extend(matchers);

/**
 * This is used if you want to use your own cookie, specifically for local logged in tests.
 * Only used if USE_COOKIE is set to true and will overwrite all other cookies that may be set.
 * Replace `value` string with your own auth cookie.
 * NOTE:: the string typically starts with "ey" and ends with "="
 */
const MANUAL_COOKIE = [
  {
    //DO NOT commit your cookie value into repo.
    value: "ey...=",
    name: "cxguser",
    domain: "api.cellxgene.dev.single-cell.czi.technology",
    path: "/",
    expires: 1965216850,
    httpOnly: true,
    secure: true,
    sameSite: "None",
  } as never,
];

/**
 * (thuang): Add `czi-checker`, so Plausible will ignore it.
 * NOTE: This changes all browsers to use Desktop Chrome UA, so please look
 * out for bugs that could be caused by this.
 * https://github.com/matomo-org/device-detector/blob/master/regexes/bots.yml#L2762
 */
const CZI_CHECKER = " czi-checker";

/**
 * Set this environment variable to enable retry
 */
const SHOULD_RETRY = process.env.RETRY !== "false";

const CLIPBOARD_PERMISSIONS = ["clipboard-read", "clipboard-write"];
if (!SHOULD_RETRY) {
  console.log('Skipping retry because "RETRY" is set to false');
}

// 'github' for GitHub Actions CI to generate annotations, default otherwise
const PLAYWRIGHT_REPORTER = process.env.CI
  ? ([["github"], ["line"], ["allure-playwright"]] as ReporterDescription[])
  : ([
      ["list"],
      [
        "html",
        {
          open: "failure",
          host: "localhost",
          port: 9220,
          outputFolder: "./html-reports",
        },
      ],
    ] as ReporterDescription[]);

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

  /* Folder for test artifacts such as screenshots, videos, traces, etc. */
  outputDir: "playwright-report/",

  /* Configure projects for major browsers */
  projects: [
    {
      name: "setup",
      testMatch: "**/*.setup.ts",
    },
    {
      name: "chromium",
      dependencies: ["setup"],
      use: {
        ...devices["Desktop Chrome"],
        userAgent: devices["Desktop Chrome"].userAgent + CZI_CHECKER,
        permissions: CLIPBOARD_PERMISSIONS,
      },
    },
    {
      name: "firefox",
      dependencies: ["setup"],
      use: {
        ...devices["Desktop Firefox"],
        userAgent: devices["Desktop Firefox"].userAgent + CZI_CHECKER,
      },
    },
    {
      name: "edge",
      dependencies: ["setup"],
      use: {
        ...devices["Desktop Edge"],
        userAgent: devices["Desktop Edge"].userAgent + CZI_CHECKER,
        permissions: CLIPBOARD_PERMISSIONS,
      },
    },
  ],

  /* Reporter to use. See https://playwright.dev/docs/test-reporters */
  reporter: PLAYWRIGHT_REPORTER,

  retries: SHOULD_RETRY ? 2 : 0,

  /* The base directory, relative to the config file, for snapshot files created with toMatchSnapshot and toHaveScreenshot. */
  snapshotDir: "./__snapshots__",

  /**
   * (thuang): For colocation, component test files live in the same directory as
   * their test pages. In the future, when Playwright component tests exit Beta,
   * we can move test files to colocate with their component files.
   *
   * https://playwright.dev/docs/test-components
   */
  // testDir: "tests",

  /**
   * Maximum time one test can run for.
   * (thuang): 5 mins because FF and Edge need extra time to tear down context
   */
  timeout: 5 * 60 * 1000,

  use: {
    ...COMMON_PLAYWRIGHT_CONTEXT,
    storageState: getStorageState(),
  },

  /* Opt out of parallel tests. */
  /**
   * By default Github Action's hosted runner has 2 CPUs, so Playwright will
   * spin up NUM_CPU / 2 workers, which is 1 worker. But locally it will run
   * tests with more workers
   * https://github.com/microsoft/playwright/issues/19408#issuecomment-1347341819
   */
  // workers: 1,

  /* Run your local dev server before starting the tests */
  // webServer: {
  //   command: 'npm run start',
  //   port: 3000,
  // },
  maxFailures: process.env.CI ? CICD_MAX_FAILURE : undefined,
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
  if (fs.existsSync(LOGIN_STATE_FILENAME)) {
    const loginState = JSON.parse(
      fs.readFileSync(LOGIN_STATE_FILENAME, "utf-8")
    );

    const storageState = getFeatureFlags({ curator: true });
    // Merge loginState with featureFlags
    storageState.cookies = [...storageState.cookies, ...loginState.cookies];
    storageState.origins = [...storageState.origins, ...loginState.origins];

    return storageState;
  }

  // For testing auth tests locally with a manual cookie
  if (process.env.USE_COOKIE === "true") {
    const storageState = getFeatureFlags({ curator: true });
    storageState.cookies = MANUAL_COOKIE;

    return storageState;
  }

  return getFeatureFlags();
}

export default defineConfig(config);
