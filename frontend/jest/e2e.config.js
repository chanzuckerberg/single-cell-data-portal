const nextJest = require("next/jest");

const createJestConfig = nextJest({
  // Provide the path to your Next.js app to load next.config.js and .env files in your test environment
  dir: "./",
});

module.exports = createJestConfig({
  globalSetup: "./jest/playwright-globalSetup.js",
  moduleDirectories: ["node_modules", "<rootDir>"],
  moduleNameMapper: {
    "\\.svg$": "<rootDir>/jest/__mocks__/svgMock.js",
  },
  preset: "jest-playwright-preset",
  rootDir: "../",
  setupFilesAfterEnv: ["expect-playwright", "./jest/playwright.setup.js"],
  testMatch: ["**/tests/**/?(*.)(spec|test).[jt]s?(x)"],
  testRunner: "jest-circus/runner",
});
