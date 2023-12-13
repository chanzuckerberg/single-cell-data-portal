import {
  CommonThemeProps,
  getColors,
  getCorners,
  getFontWeights,
  getPalette,
  getShadows,
  getSpaces,
} from "@czi-sds/components";
import { createTheme } from "@mui/material/styles";
import { defaultAppTheme, makeThemeOptions } from "@czi-sds/components";

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

export const spacesXxl = (props: CommonThemeProps) => getSpaces(props)?.xxl;
export const spacesXl = (props: CommonThemeProps) => getSpaces(props)?.xl;
export const spacesL = (props: CommonThemeProps) => getSpaces(props)?.l;
export const spacesM = (props: CommonThemeProps) => getSpaces(props)?.m;
export const spacesS = (props: CommonThemeProps) => getSpaces(props)?.s;
export const spacesXs = (props: CommonThemeProps) => getSpaces(props)?.xs;
export const spacesXxs = (props: CommonThemeProps) => getSpaces(props)?.xxs;
export const spacesXxxs = (props: CommonThemeProps) => getSpaces(props)?.xxxs;
export const spacesDefault = (props: CommonThemeProps) =>
  getSpaces(props)?.default;

const corners = {
  l: 16,
  m: 4,
  none: 0,
  s: 2,
};

export const cornersL = (props: CommonThemeProps) => getCorners(props)?.l;
export const cornersM = (props: CommonThemeProps) => getCorners(props)?.m;
export const cornersS = (props: CommonThemeProps) => getCorners(props)?.s;
export const cornersNone = (props: CommonThemeProps) => getCorners(props)?.none;

const typography = {
  fontFamily: "Inter",
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
        letterSpacing: "-0.08px",
        lineHeight: "20px",
      },
      xs: {
        fontSize: 13,
        fontWeight: fontWeights.regular,
        letterSpacing: "-0.04px",
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
        letterSpacing: "-0.05px",
        lineHeight: "16px",
      },
    },
    caps: {
      xxs: {
        fontSize: 12,
        fontWeight: fontWeights.semibold,
        letterSpacing: "0.36px",
        lineHeight: "16px",
        textTransform: "uppercase" as const,
      },
      xxxs: {
        fontSize: 11,
        fontWeight: fontWeights.semibold,
        letterSpacing: "0.33px",
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
        fontSize: 18,
        fontWeight: fontWeights.semibold,
        letterSpacing: "-0.31px",
        lineHeight: "24px",
      },
      m: {
        fontSize: 16,
        fontWeight: fontWeights.semibold,
        letterSpacing: "-0.18px",
        lineHeight: "20px",
      },
      s: {
        fontSize: 14,
        fontWeight: fontWeights.semibold,
        letterSpacing: "-0.1px",
        lineHeight: "20px",
      },
      xl: {
        fontSize: 24,
        fontWeight: fontWeights.semibold,
        letterSpacing: "-0.37px",
        lineHeight: "32px",
      },
      xs: {
        fontSize: 13,
        fontWeight: fontWeights.semibold,
        letterSpacing: "0px",
        lineHeight: "16px",
      },
      xxl: {
        fontSize: 32,
        fontWeight: fontWeights.semibold,
        letterSpacing: "-0.56px",
        lineHeight: "36px",
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

// Colors

const primaryColors = {
  "100": "#EBF5FF",
  "200": "#7DBCFF",
  "300": "#4599FF",
  "400": "#0073FF",
  "500": "#0056C6",
  "600": "#00429F",
};

export const textPrimary = (props: CommonThemeProps) =>
  getPalette(props)?.text?.primary;

export const textSecondary = (props: CommonThemeProps) =>
  getPalette(props)?.text?.secondary;

export const primary100 = (props: CommonThemeProps) =>
  getColors(props)?.primary[100];

export const primary400 = (props: CommonThemeProps) =>
  getColors(props)?.primary[400];

export const primary500 = (props: CommonThemeProps) =>
  getColors(props)?.primary[500];

export const primary600 = (props: CommonThemeProps) =>
  getColors(props)?.primary[600];

const infoColors = {
  "100": primaryColors["100"],
  "200": primaryColors["200"],
  "400": primaryColors["400"],
  "600": primaryColors["600"],
};

themeOptions.colors.primary = primaryColors;
themeOptions.colors.info = infoColors;

export const success100 = (props: CommonThemeProps) =>
  getColors(props)?.success[100];

export const success500 = (props: CommonThemeProps) =>
  getColors(props)?.success[500];

export const warning100 = (props: CommonThemeProps) =>
  getColors(props)?.warning[100];

export const warning400 = (props: CommonThemeProps) =>
  getColors(props)?.warning[400];

export const warning500 = (props: CommonThemeProps) =>
  getColors(props)?.warning[500];

export const error100 = (props: CommonThemeProps) =>
  getColors(props)?.error[100];

export const error400 = (props: CommonThemeProps) =>
  getColors(props)?.error[400];

export const error500 = (props: CommonThemeProps) =>
  getColors(props)?.error[500];

export const grey100 = (props: CommonThemeProps) => getColors(props)?.gray[100];
export const gray100 = grey100;

export const grey200 = (props: CommonThemeProps) => getColors(props)?.gray[200];
export const gray200 = grey200;

export const grey300 = (props: CommonThemeProps) => getColors(props)?.gray[300];
export const gray300 = grey300;

export const grey400 = (props: CommonThemeProps) => getColors(props)?.gray[400];
export const gray400 = grey400;

export const grey500 = (props: CommonThemeProps) => getColors(props)?.gray[500];
export const gray500 = grey500;

export const grey600 = (props: CommonThemeProps) => getColors(props)?.gray[600];
export const gray600 = grey600;

export const greyWhite = () => "#ffffff";
export const grayWhite = greyWhite;

themeOptions.colors.gray = { ...themeOptions.colors.gray, "400": "#999999" };

export const beta100 = (props: CommonThemeProps) => getColors(props)?.beta[100];

export const beta400 = (props: CommonThemeProps) => getColors(props)?.beta[400];

export const beta600 = (props: CommonThemeProps) => getColors(props)?.beta[600];

export const OFF_WHITE = "#f8f8f8";

export const PINK = "#E9429A";

// Font Weights
export const fontWeightBold = (props: CommonThemeProps) =>
  getFontWeights(props)?.bold;

export const fontWeightLight = (props: CommonThemeProps) =>
  getFontWeights(props)?.light;

export const fontWeightMedium = (props: CommonThemeProps) =>
  getFontWeights(props)?.medium;

export const fontWeightRegular = (props: CommonThemeProps) =>
  getFontWeights(props)?.regular;

/**
 * font-weight 600
 */
export const fontWeightSemibold = (props: CommonThemeProps) =>
  getFontWeights(props)?.semibold;

// Shadow
export const shadowL = (props: CommonThemeProps) => getShadows(props)?.l;
export const shadowM = (props: CommonThemeProps) => getShadows(props)?.m;
export const shadowS = (props: CommonThemeProps) => getShadows(props)?.s;

const appTheme = makeThemeOptions(themeOptions);

export const theme = createTheme(appTheme);
