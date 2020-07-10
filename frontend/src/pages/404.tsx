import React, { FC } from "react";
import Layout from "../components/layout";
import SEO from "../components/seo";

const NotFoundPage: FC = () => (
  <Layout>
    <SEO title="404: Not found" />
    <h1>NOT FOUND</h1>
    <p>We are still in the lab sequencing that one... </p>
  </Layout>
);

export default NotFoundPage;
