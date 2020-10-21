import React, { FC } from "react";
import { isSSR } from "src/common/utils/isSSR";
import AppContainer from "src/components/AppContainer";
import CookieBanner from "src/components/CookieBanner";
import Layout from "../components/Layout";
import SEO from "../components/seo";

const Index: FC = () => {
  if (isSSR()) return null;

  return (
    <Layout>
      <SEO title="Explore Data" />
      <AppContainer />
      <CookieBanner />
    </Layout>
  );
};

export default Index;
