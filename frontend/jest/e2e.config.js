module.exports = {
  moduleDirectories: ["node_modules", "."],
  preset: "jest-playwright-preset",
  rootDir: "../",
  setupFilesAfterEnv: ["expect-playwright", "./jest/playwright.setup.js"],
  testMatch: ["**/tests/**/?(*.)(spec|test).[jt]s?(x)"],
  testRunner: "jest-circus/runner",
  transform: {
    "^.+\\.(tsx?|jsx?)$": `<rootDir>/jest/jest-preprocess.js`,
  },
};
