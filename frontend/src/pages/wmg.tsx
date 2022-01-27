import { ThemeProvider as EmotionThemeProvider } from "@emotion/react";
import { StylesProvider, ThemeProvider } from "@material-ui/core/styles";
import createTheme from "@material-ui/core/styles/createTheme";
import { defaultAppTheme, makeThemeOptions } from "czifui";
import WheresMyGene from "src/views/WheresMyGene";

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

const Page = (): JSX.Element => (
  <StylesProvider injectFirst>
    <EmotionThemeProvider theme={theme}>
      <ThemeProvider theme={theme}>
        <WheresMyGene />
      </ThemeProvider>
    </EmotionThemeProvider>
  </StylesProvider>
);

export default Page;
