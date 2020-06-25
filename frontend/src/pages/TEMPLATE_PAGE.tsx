import React, { FC } from "react";
import { Link } from "gatsby";

import Layout from "../components/layout";
import SEO from "../components/seo";

const Template: FC = () => (
  <Layout>
    <SEO title="Template" />
    <h1>Foo</h1>
    <p>Bar</p>
    <Link to="/">url to home</Link>
  </Layout>
);

export default Template;
