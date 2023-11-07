import { ThemeProvider as EmotionThemeProvider } from "@emotion/react";
import { StyledEngineProvider, ThemeProvider } from "@mui/material/styles";
import { datadogRum } from "@datadog/browser-rum";
import { NextPage } from "next";
import { AppProps } from "next/app";
import Head from "next/head";
import Script from "next/script";
import { FC, useEffect } from "react";
import { QueryClient, QueryClientProvider } from "react-query";
import { ReactQueryDevtools } from "react-query/devtools";
import { EVENTS } from "src/common/analytics/events";
import { checkFeatureFlags } from "src/common/featureFlags";
import { networkGuard } from "src/common/networkGuard";
import { theme } from "src/common/theme";
import DefaultLayout from "src/components/Layout/components/defaultLayout";
import configs from "src/configs/configs";
import "src/global.scss";
// (thuang): `layout.css` needs to be imported after `global.scss`
import "src/layout.css";

const OG_PAGE_TITLE = "Cellxgene Data Portal";

declare global {
  interface Window {
    plausible: {
      q: unknown[];
      (event: EVENTS, options?: { props: { [key: string]: unknown } }): void;
    };
  }
}

const queryClient = new QueryClient();

checkFeatureFlags();
networkGuard();

type NextPageWithLayout = NextPage & {
  Layout?: FC;
};

type AppPropsWithLayout = AppProps & {
  Component: NextPageWithLayout;
};

datadogRum.init({
  applicationId: "44f77ca2-1482-404a-ad38-23499bb925e5",
  clientToken: "pub55d4baaac2091f9656a83da732732a89",
  site: "datadoghq.com",
  service: "single-cell-data-portal",
  env: process.env.NODE_ENV,
  sessionSampleRate: 100,
  sessionReplaySampleRate: 20,
  trackUserInteractions: true,
  trackResources: true,
  trackLongTasks: true,
  defaultPrivacyLevel: "mask-user-input",
  allowedTracingUrls: [
    "<https://pr-6179-backend.rdev.single-cell.czi.technology>",
    /https:\/\/.*\.my-api-domain\.com/,
    (url) =>
      url.startsWith(
        "<https://pr-6179-backend.rdev.single-cell.czi.technology>"
      ),
  ],
});

function App({ Component, pageProps }: AppPropsWithLayout): JSX.Element {
  const Layout = Component.Layout || DefaultLayout;

  // (thuang): Per Plausible instruction
  // "...make sure your tracking setup includes the second line as shown below"
  // https://plausible.io/docs/custom-event-goals#1-trigger-custom-events-with-javascript-on-your-site
  useEffect(() => {
    window.plausible = window.plausible || tempPlausible;

    function tempPlausible(...args: unknown[]): void {
      (window.plausible.q = window.plausible.q || []).push(args);
    }
  }, []);

  return (
    <>
      <Head>
        {/* meta tag must have a `key` prop to allow page specific overwrites */}
        {/* See: https://nextjs.org/docs/pages/api-reference/components/head */}
        <meta name="twitter:card" key="twitter:card" content="summary" />

        {/* Open Graph */}
        <meta
          property="og:url"
          key="og:url"
          content="https://cellxgene.cziscience.com/"
        />
        <meta
          property="og:image"
          key="og:image"
          content={"https://cellxgene.cziscience.com/open-graph.jpg"}
        />
        <meta property="og:creator" key="og:creator" content="@cziscience" />
        <meta property="og:site" key="og:site" content="@cziscience" />
        <meta
          property="og:site_name"
          key="og:site_name"
          content={OG_PAGE_TITLE}
        />
        <meta property="og:title" key="og:title" content={OG_PAGE_TITLE} />
        <meta
          property="og:description"
          key="og:description"
          content="Find, download, and visually explore curated and standardized single cell datasets."
        />
        <meta
          property="og:image"
          key="og:image"
          content="https://cellxgene.cziscience.com/open-graph.jpg"
        />
        <meta
          property="twitter:image"
          key="twitter:image"
          content="https://cellxgene.cziscience.com/open-graph.jpg"
        />

        <meta property="twitter:card" key="twitter:card" content="summary" />
        <meta property="og:type" key="og:type" content="website" />
      </Head>
      <QueryClientProvider client={queryClient}>
        <StyledEngineProvider>
          <EmotionThemeProvider theme={theme}>
            <ThemeProvider theme={theme}>
              <Layout>
                <Component {...pageProps} />
              </Layout>
              <ReactQueryDevtools />
            </ThemeProvider>
          </EmotionThemeProvider>
        </StyledEngineProvider>
      </QueryClientProvider>
      <Script
        data-domain={configs.PLAUSIBLE_DATA_DOMAIN}
        src="https://plausible.io/js/script.js"
      />
    </>
  );
}

export default App;
