import { NextPage } from "next";
import { AppProps } from "next/app";
import Script from "next/script";
import { FC } from "react";
import { QueryClient, QueryClientProvider } from "react-query";
import { ReactQueryDevtools } from "react-query/devtools";
import { checkFeatureFlags } from "src/common/featureFlags";
import { FEATURES } from "src/common/featureFlags/features";
import { useFeatureFlag } from "src/common/hooks/useFeatureFlag";
import DefaultLayout from "src/components/Layout/components/defaultLayout";
import configs from "src/configs/configs";
import "src/global.scss";
// (thuang): `layout.css` needs to be imported after `global.scss`
import "src/layout.css";

const queryClient = new QueryClient();

checkFeatureFlags();

type NextPageWithLayout = NextPage & {
  Layout?: FC;
};

type AppPropsWithLayout = AppProps & {
  Component: NextPageWithLayout;
};

function App({ Component, pageProps }: AppPropsWithLayout): JSX.Element {
  const isFilterEnabled = useFeatureFlag(FEATURES.FILTER);
  const Layout = (isFilterEnabled && Component.Layout) || DefaultLayout;
  return (
    <>
      <QueryClientProvider client={queryClient}>
        <Layout>
          <Component {...pageProps} />
        </Layout>
        <ReactQueryDevtools />
      </QueryClientProvider>
      <Script
        data-domain={configs.PLAUSIBLE_DATA_DOMAIN}
        src="https://plausible.io/js/plausible.js"
      />
    </>
  );
}

export default App;
