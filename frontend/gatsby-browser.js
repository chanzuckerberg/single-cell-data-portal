/**
 * Implement Gatsby's Browser APIs in this file.
 *
 * See: https://www.gatsbyjs.org/docs/browser-apis/
 */

import React from "react";
import { Auth0Provider } from "./src/contexts/auth0Context";

export const wrapRootElement = ({ element }) => {
  return <Auth0Provider>{element}</Auth0Provider>;
};
