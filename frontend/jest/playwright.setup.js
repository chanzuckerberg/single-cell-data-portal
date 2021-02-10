/**
 * `frontend/jest-playwright.config.js` is for configuring Playwright's launch config options
 * `frontend/jest/playwright.setup.js` is for configuring `jest`, `browser`,
 * and `page` objects
 */

// (thuang): This is the max time a test can take to run.
// Since when debugging, we run and !headless, this means
// a test can take more time to finish, so we don't want
// jest to shut off the test too soon
jest.setTimeout(2 * 60 * 1000);

// (thuang): Please make sure this number matches
// `RETRY_ATTEMPTS` in `jest/screenshot_env.js` (when we add screenshot)
jest.retryTimes(2);

beforeEach(async () => {
  const client = await page.context().newCDPSession(page);

  await client.send("Animation.setPlaybackRate", {
    // (thuang): Max speed to "disable" animation
    playbackRate: 12,
  });
});
