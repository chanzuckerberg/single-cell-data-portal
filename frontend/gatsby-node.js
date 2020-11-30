/**
 * Implement Gatsby's Node APIs in this file.
 *
 * See: https://www.gatsbyjs.org/docs/node-apis/
 */
exports.onCreatePage = ({ page, actions }) => {
  const { createPage } = actions;

  // Make the index page match everything client side.
  if (page.path === `/`) {
    page.matchPath = `/*`;
    createPage(page);
  }
};
