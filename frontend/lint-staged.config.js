module.exports = {
  "*.js": "prettier --write",
  "*.ts?(x)": filenames =>
    `concurrently "prettier --parser typescript --write ${filenames.join(
      " "
    )}" "npm run check-types"`,
};
