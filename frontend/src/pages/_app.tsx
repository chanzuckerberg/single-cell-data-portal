import { Auth0Provider } from '@auth0/auth0-react';
import { ThemeProvider as EmotionThemeProvider } from "@emotion/react";
import { StyledEngineProvider, ThemeProvider } from "@mui/material/styles";
import { NextPage } from "next";
import { AppProps } from "next/app";
import Script from "next/script";
import { FC, useEffect, useState } from "react";
import { QueryClient, QueryClientProvider } from "react-query";
import { ReactQueryDevtools } from "react-query/devtools";
import { EVENTS } from "src/common/analytics/events";
import { checkFeatureFlags } from "src/common/featureFlags";
import { theme } from "src/common/theme";
import DefaultLayout from "src/components/Layout/components/defaultLayout";
import configs from "src/configs/configs";
import "src/global.scss";
// (thuang): `layout.css` needs to be imported after `global.scss`
import "src/layout.css";

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

type NextPageWithLayout = NextPage & {
  Layout?: FC;
};

type AppPropsWithLayout = AppProps & {
  Component: NextPageWithLayout;
};

function App({ Component, pageProps }: AppPropsWithLayout): JSX.Element {
  const Layout = Component.Layout || DefaultLayout;

  // (thuang): Per Plausible instruction
  // "...make sure your tracking setup includes the second line as shown below"
  // https://plausible.io/docs/custom-event-goals#1-trigger-custom-events-with-javascript-on-your-site
  useEffect(() => {
    console.log("PLAUSIBLE");
    window.plausible = window.plausible || tempPlausible;

    function tempPlausible(...args: unknown[]): void {
      (window.plausible.q = window.plausible.q || []).push(args);
    }
  }, []);
  console.log("main body");
  let windowLocationOrigin;

  const [isMounted, setMounted] = useState(false);

  useEffect(() => {
    console.log("RERENDERED");
    console.log(`query client ${queryClient}`); 
    windowLocationOrigin = window.location.origin;
    setMounted(true);
  }, []);
  console.log(`query client ${queryClient}`);
  return (
    isMounted ? (
      <>
        <Auth0Provider
          domain="corporanet.local"
          clientId="local-client-id"
          redirectUri={windowLocationOrigin}
        >
          <QueryClientProvider client={queryClient}>
            <StylesProvider injectFirst>
              <EmotionThemeProvider theme={theme}>
                <ThemeProvider theme={theme}>
                  <Layout>
                    <Component {...pageProps} />
                  </Layout>
                  <ReactQueryDevtools />
                </ThemeProvider>
              </EmotionThemeProvider>
            </StylesProvider>
          </QueryClientProvider>
        </Auth0Provider>
        <Script
          data-domain={configs.PLAUSIBLE_DATA_DOMAIN}
          src="https://plausible.io/js/script.js"
        />
      </>
    ) : (<></>)
  );
}

export default App;
