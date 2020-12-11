import React, { FC } from "react";
import { QueryCache, ReactQueryCacheProvider } from "react-query";
import { ReactQueryDevtools } from "react-query-devtools";
import { isSSR } from "src/common/utils/isSSR";
import AppContainer from "src/components/AppContainer";
import CookieBanner from "src/components/CookieBanner";
import Layout from "../components/Layout";
import SEO from "../components/seo";

const queryCache = new QueryCache();

const Index: FC = () => {
  if (isSSR()) return null;

  return (
    <ReactQueryCacheProvider queryCache={queryCache}>
      <Layout>
        <SEO title="Explore Data" />
        <AppContainer />
        <CookieBanner />
      </Layout>
      <ReactQueryDevtools />
    </ReactQueryCacheProvider>
  );
};

export default Index;
