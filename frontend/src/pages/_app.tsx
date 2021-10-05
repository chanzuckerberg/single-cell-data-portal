import { AppProps } from "next/app";
import Script from "next/script";
import { QueryCache, ReactQueryCacheProvider } from "react-query";
import { ReactQueryDevtools } from "react-query-devtools";
import { checkFeatureFlags } from "src/common/featureFlags";
import configs from "src/configs/configs";
import "src/global.scss";
// (thuang): `layout.css` needs to be imported after `global.scss`
import "src/layout.css";
import Layout from "../components/Layout";

const queryCache = new QueryCache();

checkFeatureFlags();

function App({ Component, pageProps }: AppProps): JSX.Element {
  return (
    <>
      <ReactQueryCacheProvider queryCache={queryCache}>
        <Layout>
          <Component {...pageProps} />
        </Layout>
        <ReactQueryDevtools />
      </ReactQueryCacheProvider>
      <Script
        data-domain={configs.PLAUSIBLE_DATA_DOMAIN}
        src="https://plausible.io/js/plausible.js"
      />
    </>
  );
}

export default App;
