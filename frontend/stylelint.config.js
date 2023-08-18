module.exports = {
  extends: "stylelint-config-recommended",
  ignoreFiles: [
    "public/**/*",
    // (thuang): Ignore `venv` folder
    "venv/**/*",
    "html-reports/**/*",
  ],
};
