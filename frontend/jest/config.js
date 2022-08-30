const nextJest = require("next/jest");

const createJestConfig = nextJest({
  // Provide the path to your Next.js app to load next.config.js and .env files in your test environment
  dir: "./",
});

module.exports = createJestConfig({
  coverageDirectory: "<rootDir>/client-coverage",
  coveragePathIgnorePatterns: ["<rootDir>/node_modules/", "<rootDir>/build/"],
  coverageReporters: ["text-summary", "json", "html"],
  coverageThreshold: {
    global: {
      branches: 35,
      functions: 40,
      lines: 55,
      statements: 55,
    },
  },
  globals: {},
  moduleDirectories: ["node_modules", "src"],
  moduleFileExtensions: ["ts", "tsx", "js", "jsx"],
  moduleNameMapper: {},
  modulePaths: ["<rootDir>/"],
  rootDir: "../",
  setupFiles: ["<rootDir>/jest/test-setup.js"],
  testPathIgnorePatterns: [
    "<rootDir>/node_modules/",
    "<rootDir>/build/",
    "<rootDir>/tasks/",
    "<rootDir>/src/server",
    "<rootDir>/cypress/",
    "<rootDir>/puppeteer/",
  ],
  testRegex: "(/__tests__/.*|(\\.|/)(test|spec))\\.(tsx?|jsx?)$",
  transform: {
    "^.+\\.(js|ts)x?$": "babel-jest",
  },
});
