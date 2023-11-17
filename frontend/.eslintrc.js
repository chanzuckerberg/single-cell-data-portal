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
    "plugin:@typescript-eslint/recommended",
    "prettier",
    "plugin:prettier/recommended",
    "plugin:sonarjs/recommended",
    "plugin:mdx/recommended",
    "next",
  ],
  overrides: [
    // Override some TypeScript rules just for .js files
    {
      files: ["*.js"],
      rules: {
        "@typescript-eslint/no-var-requires": "off",
      },
    },
  ],
  parser: "@typescript-eslint/parser",
  parserOptions: {
    project: "./tsconfig.json",
    ecmaFeatures: {
      jsx: true,
    },
    ecmaVersion: 2020,
    // Allows for the parsing of modern ECMAScript features
    sourceType: "module", // Allows for the use of imports
  },
  plugins: ["@typescript-eslint", "sonarjs", "@blueprintjs", "jsx-expressions"],
  rules: {
    "@typescript-eslint/camelcase": 0,
    // Disable prop-types as we use TypeScript for type checking
    "@typescript-eslint/explicit-function-return-type": "off",
    // (thuang): Allow args prefixed with `_`
    // example: https://eslint.org/docs/rules/no-unused-vars#argsignorepattern
    "@typescript-eslint/no-unused-vars": [
      "error",
      {
        args: "after-used",
        argsIgnorePattern: "^_",
        ignoreRestSiblings: false,
        vars: "all",
      },
    ],
    camelcase: "off",
    "react/jsx-no-target-blank": 0,
    "react/prop-types": "off",
    // (thuang): We use nested template literals extensively
    "sonarjs/no-nested-template-literals": "off",
    // React Hooks
    "react-hooks/rules-of-hooks": "error",
    "react-hooks/exhaustive-deps": "error",
    "jsx-expressions/strict-logical-expressions": "error",
  },
  settings: {
    "mdx/code-blocks": true,
    react: {
      version: "detect",
    },
  },
};
