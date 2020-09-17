/**
 * `frontend/jest-playwright.config.js` is for configuring Playwright's launch config options
 * `frontend/jest/playwright.setup.js` is for configuring `jest`, `browser`,
 * and `page` objects
 */

const isHeadful =
  process.env.HEADFUL === "true" || process.env.HEADLESS === "false";

const DEFAULT_LAUNCH_CONFIG = {
  args: ["--ignore-certificate-errors", "--ignore-ssl-errors"],
  headless: !isHeadful,
  ignoreHTTPSErrors: true,
};

const DEFAULT_CONTEXT_CONFIG = {
  acceptDownloads: true,
};

module.exports = {
  browserContext: "incognito",
  contextOptions: DEFAULT_CONTEXT_CONFIG,
  launchOptions: DEFAULT_LAUNCH_CONFIG,
};
