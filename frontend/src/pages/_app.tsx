import { ThemeProvider as EmotionThemeProvider } from "@emotion/react";
import { StylesProvider, ThemeProvider } from "@material-ui/core/styles";
import createTheme from "@material-ui/core/styles/createTheme";
import { defaultAppTheme, makeThemeOptions } from "czifui";
import { NextPage } from "next";
import { AppProps } from "next/app";
import Script from "next/script";
import { FC } from "react";
import { QueryClient, QueryClientProvider } from "react-query";
import { ReactQueryDevtools } from "react-query/devtools";
import { checkFeatureFlags } from "src/common/featureFlags";
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

const customTheme = {
  typography: {
    ...defaultAppTheme.typography,
    fontFamily: "Roboto",
  },
};

const themeOptions = { ...defaultAppTheme, ...customTheme };

const primaryColors = {
  "100": "#dfefff",
  "200": "#73b7ff",
  "300": "#459fff",
  "400": "#0e7dec",
  "500": "#005fc6",
  "600": "#004c9f",
};

const infoColors = {
  "100": "#73b7ff",
  "200": "#F3EDFC",
  "400": "#0e7dec",
  "600": "#004c9f",
};

themeOptions.colors.primary = primaryColors;
themeOptions.colors.info = infoColors;

const appTheme = makeThemeOptions(themeOptions);

const theme = createTheme(appTheme);

function App({ Component, pageProps }: AppPropsWithLayout): JSX.Element {
  const Layout = Component.Layout || DefaultLayout;
  return (
    <>
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
      <Script
        data-domain={configs.PLAUSIBLE_DATA_DOMAIN}
        src="https://plausible.io/js/plausible.js"
      />
    </>
  );
}

export default App;
