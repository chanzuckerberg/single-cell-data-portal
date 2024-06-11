import {
  CommonThemeProps,
  SDSAppTheme,
  getColors,
  getCorners,
  getFontWeights,
  getPalette,
  getShadows,
  getSpaces,
  makeThemeOptions,
} from "@czi-sds/components";
import { createTheme } from "@mui/material/styles";

import { Inter, IBM_Plex_Mono } from "next/font/google";

export const INTER_FONT_CSS_VARIABLE = "--font-inter";

/**
 * (thuang): We should only load the font once in the app.
 * If loading a variable font, you don't need to specify the font weight
 * `Inter` is a variable font
 */
export const inter = Inter({
  subsets: ["latin"],
  /**
   * (thuang): We can't use the font-weight variable here, according to Next.js warning
   * CSS variable here is used in CSS files
   */
  variable: "--font-inter",
});

/**
 * (masoudmanson): IMB_Plex_Mono is needed for SDS code blocks.
 */
export const ibm_plex_mono = IBM_Plex_Mono({
  subsets: ["latin"],
  variable: "--font-ibm-plex-mono",
  weight: ["400", "600"],
});

const { fontWeights } = SDSAppTheme;

const iconSizes = {
  input: { height: 16, width: 16 }, // for use with input icons only (radio and checkbox)
  l: { height: 24, width: 24 },
  m: { height: 16, width: 16 }, // (masoudmanson): SDS doesn't have a medium icon size!
  s: { height: 16, width: 16 }, // (masoudmanson): Used to be 12px
  xl: { height: 32, width: 32 },
  xs: { height: 12, width: 12 }, // (masoudmanson): Used to be 8px
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

/**
 * (masoudmanson): SDS has introduced new font styles for tabular numbers.
 */
const tabularNums = "tabular-nums";

const typography = {
  fontFamily: {
    body: inter.style.fontFamily,
    caps: inter.style.fontFamily,
    code: ibm_plex_mono.style.fontFamily,
    header: inter.style.fontFamily,
    tabular: inter.style.fontFamily,
  },
  styles: {
    body: {
      regular: {
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
      semibold: {
        button: {
          fontSize: 14,
          fontWeight: fontWeights.semibold,
          letterSpacing: "0px",
          lineHeight: "20px",
          textTransform: "none" as const,
        },
        l: {
          fontSize: 18,
          fontWeight: fontWeights.semibold,
          letterSpacing: "0px",
          lineHeight: "24px",
        },
        m: {
          fontSize: 16,
          fontWeight: fontWeights.semibold,
          letterSpacing: "0px",
          lineHeight: "24px",
        },
        s: {
          fontSize: 14,
          fontWeight: fontWeights.semibold,
          letterSpacing: "-0.08px",
          lineHeight: "20px",
        },
        xs: {
          fontSize: 13,
          fontWeight: fontWeights.semibold,
          letterSpacing: "-0.04px",
          lineHeight: "20px",
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
          letterSpacing: "-0.05px",
          lineHeight: "16px",
        },
      },
    },
    caps: {
      semibold: {
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
    },
    code: {
      regular: {
        s: {
          fontSize: 14,
          fontWeight: fontWeights.regular,
          letterSpacing: "0px",
          lineHeight: "24px",
          textTransform: "none" as const,
        },
        xs: {
          fontSize: 13,
          fontWeight: fontWeights.regular,
          letterSpacing: "0px",
          lineHeight: "20px",
          textTransform: "none" as const,
        },
      },
      semibold: {
        s: {
          fontSize: 14,
          fontWeight: fontWeights.semibold,
          letterSpacing: "0px",
          lineHeight: "24px",
          textTransform: "none" as const,
        },
        xs: {
          fontSize: 13,
          fontWeight: fontWeights.semibold,
          letterSpacing: "0px",
          lineHeight: "20px",
          textTransform: "none" as const,
        },
      },
    },
    header: {
      semibold: {
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
    tabular: {
      regular: {
        s: {
          fontSize: 14,
          fontVariantNumeric: tabularNums,
          fontWeight: fontWeights.regular,
          letterSpacing: "0px",
          lineHeight: "24px",
          textTransform: "none" as const,
        },
        xs: {
          fontSize: 13,
          fontVariantNumeric: tabularNums,
          fontWeight: fontWeights.regular,
          letterSpacing: "0px",
          lineHeight: "20px",
          textTransform: "none" as const,
        },
      },
      semibold: {
        s: {
          fontSize: 14,
          fontVariantNumeric: tabularNums,
          fontWeight: fontWeights.semibold,
          letterSpacing: "0px",
          lineHeight: "24px",
          textTransform: "none" as const,
        },
        xs: {
          fontSize: 13,
          fontVariantNumeric: tabularNums,
          fontWeight: fontWeights.semibold,
          letterSpacing: "0px",
          lineHeight: "20px",
          textTransform: "none" as const,
        },
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

const themeOptions = { ...SDSAppTheme, ...customTheme };

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
  getColors(props)?.blue[100];

export const primary200 = (props: CommonThemeProps) =>
  getColors(props)?.blue[200];

export const primary300 = (props: CommonThemeProps) =>
  getColors(props)?.blue[300];

export const primary400 = (props: CommonThemeProps) =>
  getColors(props)?.blue[400];

export const primary500 = (props: CommonThemeProps) =>
  getColors(props)?.blue[500];

export const primary600 = (props: CommonThemeProps) =>
  getColors(props)?.blue[600];

themeOptions.colors.blue = primaryColors;

/**
 * (masoudmanson): SDS Theme does not include info colors as part of the primitive colors.
 * Primitive colors are: blue, green, red, yellow, purple, and gray.
 * Additionally, there are Semantic Text and Semantic Component colors available.
 *
 * You can browse the full list of SDS theme colors here:
 * https://chanzuckerberg.github.io/sci-components/?path=/story/bases-colors--primitive-colors
 */

// const infoColors = {
//   "100": primaryColors["100"],
//   "200": primaryColors["200"],
//   "400": primaryColors["400"],
//   "600": primaryColors["600"],
// };

// themeOptions.colors.info = infoColors;

export const success100 = (props: CommonThemeProps) =>
  getColors(props)?.green[100];

export const success400 = (props: CommonThemeProps) =>
  getColors(props)?.green[400];

export const success500 = (props: CommonThemeProps) =>
  getColors(props)?.green[500];

export const success600 = (props: CommonThemeProps) =>
  getColors(props)?.green[600];

export const warning100 = (props: CommonThemeProps) =>
  getColors(props)?.yellow[100];

export const warning400 = (props: CommonThemeProps) =>
  getColors(props)?.yellow[400];

export const warning500 = (props: CommonThemeProps) =>
  getColors(props)?.yellow[500];

export const error100 = (props: CommonThemeProps) => getColors(props)?.red[100];

export const error400 = (props: CommonThemeProps) => getColors(props)?.red[400];

export const error500 = (props: CommonThemeProps) => getColors(props)?.red[500];

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

export const beta100 = (props: CommonThemeProps) =>
  getColors(props)?.purple[100];

export const beta400 = (props: CommonThemeProps) =>
  getColors(props)?.purple[400];

export const beta600 = (props: CommonThemeProps) =>
  getColors(props)?.purple[600];

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
