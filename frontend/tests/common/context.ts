import { Config } from "@playwright/test";

const isHeadful =
  process.env.HEADFUL === "true" || process.env.HEADLESS === "false";

if (isHeadful) {
  console.log("Running tests in headful mode");
}

/**
 * (thuang): Keep the video size small to avoid test timing out from processing
 * large video files.
 */
const VIEWPORT = {
  height: 800,
  width: 1280,
};

export const COMMON_PLAYWRIGHT_CONTEXT: Config["use"] = {
  acceptDownloads: true,
  /* Maximum time each action such as `click()` can take. Defaults to 0 (no limit). */
  actionTimeout: 0,
  headless: !isHeadful,
  ignoreHTTPSErrors: true,
  screenshot: "only-on-failure",
  /* Collect trace when retrying the failed test. See https://playwright.dev/docs/trace-viewer */
  trace: "retain-on-failure",
  video: {
    mode: "retain-on-failure",
    size: VIEWPORT,
  },
  viewport: VIEWPORT,
};
