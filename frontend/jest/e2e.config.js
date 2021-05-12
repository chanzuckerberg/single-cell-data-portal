module.exports = {
  globalSetup: "./jest/playwright-globalSetup.js",
  moduleDirectories: ["node_modules", "<rootDir>"],
  preset: "jest-playwright-preset",
  rootDir: "../",
  setupFilesAfterEnv: ["expect-playwright", "./jest/playwright.setup.js"],
  testMatch: ["**/tests/**/?(*.)(spec|test).[jt]s?(x)"],
  testRunner: "jest-circus/runner",
  transform: {
    "^.+\\.(tsx?|jsx?)$": `<rootDir>/jest/jest-preprocess.js`,
    "^.+\\.svg$": "jest-svg-transformer",
  },
};
