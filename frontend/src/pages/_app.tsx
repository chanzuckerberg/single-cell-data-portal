import { AppProps } from "next/app";
import React from "react";
import { QueryCache, ReactQueryCacheProvider } from "react-query";
import { ReactQueryDevtools } from "react-query-devtools";
import { checkFeatureFlags } from "src/common/featureFlags";
import CookieBanner from "src/components/CookieBanner";
import "src/global.scss";
// (thuang): `layout.css` needs to be imported after `global.scss`
import "src/layout.css";
import Layout from "../components/Layout";

const queryCache = new QueryCache();

checkFeatureFlags();

function App({ Component, pageProps }: AppProps) {
  return (
    <ReactQueryCacheProvider queryCache={queryCache}>
      <Layout>
        <Component {...pageProps} />
        <CookieBanner />
      </Layout>
      <ReactQueryDevtools />
    </ReactQueryCacheProvider>
  );
}

export default App;
