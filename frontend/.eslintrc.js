// https://www.robertcooper.me/using-eslint-and-prettier-in-a-typescript-project
module.exports = {
  env: {
    browser: true,
    es6: true,
    node: true,
  },
  // Specifies the ESLint parser
  extends: [
    "eslint:recommended",
    "plugin:react/recommended",
    "plugin:@typescript-eslint/recommended",
    "prettier/@typescript-eslint",
    "plugin:prettier/recommended",
    "react-app",
    "plugin:sonarjs/recommended",
    "plugin:jest-playwright/recommended",
  ],
  overrides: [
    // Override some TypeScript rules just for .js files
    {
      files: ["*.js"],
      rules: {
        "@typescript-eslint/no-var-requires": "off", //
      },
    },
  ],
  parser: "@typescript-eslint/parser",
  parserOptions: {
    ecmaFeatures: {
      jsx: true,
    },
    ecmaVersion: 2020,
    // Allows for the parsing of modern ECMAScript features
    sourceType: "module", // Allows for the use of imports
  },
  plugins: ["@typescript-eslint", "react", "sonarjs", "@blueprintjs"],
  rules: {
    // Disable prop-types as we use TypeScript for type checking
    "@typescript-eslint/explicit-function-return-type": "off",
    "react/jsx-no-target-blank": 0,
    "react/prop-types": "off",
    "sort-keys": [
      "error",
      "asc",
      {
        caseSensitive: true,
        minKeys: 2,
        natural: false,
      },
    ],
  },
  settings: {
    react: {
      version: "detect",
    },
  },
};
