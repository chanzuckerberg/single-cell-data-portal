import { createTheme } from "@material-ui/core/styles";
import { defaultAppTheme, makeThemeOptions } from "czifui";

const { fontWeights } = defaultAppTheme;

const iconSizes = {
  input: { height: 16, width: 16 }, // for use with input icons only (radio and checkbox)
  l: { height: 24, width: 24 },
  m: { height: 16, width: 16 },
  s: { height: 12, width: 12 },
  xl: { height: 32, width: 32 },
  xs: { height: 8, width: 8 },
};

const spacing = {
  default: 8,
  l: 16,
  m: 12,
  s: 8,
  xl: 24,
  xs: 6,
  xxl: 40,
  xxs: 4,
  xxxs: 2,
};

const corners = {
  l: 16,
  m: 4,
  none: 0,
  s: 2,
};

const typography = {
  fontFamily: "Roboto Condensed",
  styles: {
    body: {
      button: {
        fontSize: 14,
        fontWeight: fontWeights.medium,
        letterSpacing: "0px",
        lineHeight: "20px",
        textTransform: "none" as const,
      },
      l: {
        fontSize: 18,
        fontWeight: fontWeights.regular,
        letterSpacing: "0px",
        lineHeight: "24px",
      },
      m: {
        fontSize: 16,
        fontWeight: fontWeights.regular,
        letterSpacing: "0px",
        lineHeight: "24px",
      },
      s: {
        fontSize: 14,
        fontWeight: fontWeights.regular,
        letterSpacing: "0px",
        lineHeight: "20px",
      },
      xs: {
        fontSize: 13,
        fontWeight: fontWeights.regular,
        letterSpacing: "0px",
        lineHeight: "20px",
      },
      xxs: {
        fontSize: 12,
        fontWeight: fontWeights.regular,
        letterSpacing: "0px",
        lineHeight: "16px",
      },
      xxxs: {
        fontSize: 11,
        fontWeight: fontWeights.regular,
        letterSpacing: "0px",
        lineHeight: "16px",
      },
    },
    caps: {
      xxs: {
        fontSize: 12,
        fontWeight: fontWeights.semibold,
        letterSpacing: "1.0px",
        lineHeight: "16px",
        textTransform: "uppercase" as const,
      },
      xxxs: {
        fontSize: 11,
        fontWeight: fontWeights.semibold,
        letterSpacing: "1.0px",
        lineHeight: "16px",
        textTransform: "uppercase" as const,
      },
      xxxxs: {
        fontSize: 10,
        fontWeight: fontWeights.semibold,
        letterSpacing: "1.0px",
        lineHeight: "12px",
        textTransform: "uppercase" as const,
      },
    },
    header: {
      l: {
        fontSize: 20,
        fontWeight: fontWeights.semibold,
        letterSpacing: "0px",
        lineHeight: "24px",
      },
      m: {
        fontSize: 16,
        fontWeight: fontWeights.semibold,
        letterSpacing: "0px",
        lineHeight: "20px",
      },
      s: {
        fontSize: 14,
        fontWeight: fontWeights.semibold,
        letterSpacing: "0px",
        lineHeight: "20px",
      },
      xl: {
        fontSize: 24,
        fontWeight: fontWeights.semibold,
        letterSpacing: "0px",
        lineHeight: "32px",
      },
      xs: {
        fontSize: 13,
        fontWeight: fontWeights.semibold,
        letterSpacing: "0px",
        lineHeight: "16px",
      },
      xxl: {
        fontSize: 28,
        fontWeight: fontWeights.semibold,
        letterSpacing: "0px",
        lineHeight: "32px",
      },
      xxs: {
        fontSize: 12,
        fontWeight: fontWeights.semibold,
        letterSpacing: "0px",
        lineHeight: "16px",
      },
      xxxs: {
        fontSize: 11,
        fontWeight: fontWeights.semibold,
        letterSpacing: "0px",
        lineHeight: "16px",
      },
    },
  },
};

const customTheme = {
  corners,
  iconSizes,
  spacing,
  typography,
};

const themeOptions = { ...defaultAppTheme, ...customTheme };

const primaryColors = {
  "100": "#EBF5FF",
  "200": "#7DBCFF",
  "300": "#4599FF",
  "400": "#0073FF",
  "500": "#0056C6",
  "600": "#00429F",
};

const infoColors = {
  "100": primaryColors["100"],
  "200": primaryColors["200"],
  "400": primaryColors["400"],
  "600": primaryColors["600"],
};

themeOptions.colors.primary = primaryColors;
themeOptions.colors.info = infoColors;

const appTheme = makeThemeOptions(themeOptions);

export const theme = createTheme(appTheme);
