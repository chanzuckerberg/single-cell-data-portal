import styled from "@emotion/styled";
import {
  CommonThemeProps,
  fontBodyS,
  getColors,
  getCorners,
  getFontWeights,
  getSpaces,
  List as SDSList,
  ListItem,
} from "czifui";
import { css } from "@emotion/react";
import { GRAY, PT_TEXT_COLOR } from "src/components/common/theme";

const cornersM = (props: CommonThemeProps) => getCorners(props)?.m;
const grey100 = (props: CommonThemeProps) => getColors(props)?.gray[100];
const grey500 = (props: CommonThemeProps) => getColors(props)?.gray[500];
const primary400 = (props: CommonThemeProps) => getColors(props)?.primary[400];
const semibold = (props: CommonThemeProps) => getFontWeights(props)?.semibold; // font-weight 600.
const spacesS = (props: CommonThemeProps) => getSpaces(props)?.s;
const spacesXs = (props: CommonThemeProps) => getSpaces(props)?.xs;
const spacesXxs = (props: CommonThemeProps) => getSpaces(props)?.xxs;
const spacesXxxs = (props: CommonThemeProps) => getSpaces(props)?.xxxs;

export const listCss = (props: CommonThemeProps) => {
  return css`
    display: grid;
    gap: ${spacesXxs(props)}px;
  `;
};

export const listSubheaderCss = (props: CommonThemeProps) => {
  return css`
    background-color: #ffffff;
    color: ${GRAY.A};
    cursor: default;
    font-size: 11px;
    font-weight: 500;
    letter-spacing: 0.03em;
    line-height: 15px;
    margin-bottom: ${spacesXxs(props)}px;
    padding-left: ${spacesS(props)}px;
    text-transform: uppercase;
  `;
};

export const listItemCss = (props: CommonThemeProps) => {
  return css`
    ${fontBodyS(props)}
    letter-spacing: -0.006em;
    margin: 0 /* overrides margin from layout.css */;

    &:before {
      display: none; /* remove list item bullet. */
    }
  `;
};

export const listItemButtonCss = (props: CommonThemeProps) => {
  return css`
    padding: ${spacesXs(props)}px ${spacesS(props)}px;

    &:hover {
      background-color: ${grey100(props)};
    }
  `;
};

export const listItemButtonSelectedCss = (props: CommonThemeProps) => {
  return css`
    &.Mui-selected {
      background-color: transparent;

      .MuiListItemText-root {
        span:first-of-type {
          font-weight: ${semibold(props)};
        }
      }
    }
  `;
};

export const listItemTextCss = (props: CommonThemeProps) => {
  return css`
    display: flex;
    gap: ${spacesS(props)}px;
    margin: 0;

    /* List item text - "primary" */

    span:first-of-type {
      flex: 1;
    }

    /* List item count - "secondary" */

    span:last-of-type {
      color: ${grey500(props)};
    }
  `;
};

export const listItemIconCss = (props: CommonThemeProps) => {
  return css`
    align-items: center;
    box-sizing: content-box;
    color: ${primary400(props)};
    display: flex;
    height: 16px;
    justify-content: center;
    margin: 0 ${spacesS(props)}px 0 0;
    min-width: 0;
    padding: ${spacesXxxs(props)}px 0;
    width: 16px;

    .MuiSvgIcon-root {
      font-size: 16px;
    }
  `;
};

export const listItemDividerCss = css`
  background-color: ${PT_TEXT_COLOR};
  border: none;
  border-bottom: 1px;
  opacity: 0.15;
  margin: 0;
`;

export const scrollbar = (props: CommonThemeProps) => {
  const colors = getColors(props);
  return css`
    &::-webkit-scrollbar {
      width: ${spacesXxs(props)}px;
    }

    &::-webkit-scrollbar-thumb {
      background-clip: content-box;
      background-color: ${colors?.gray["200"]};
      border-radius: ${cornersM(props)}px;

      &:hover {
        background-color: ${colors?.gray["300"]};
      }
    }
  `;
};

export const List = styled(SDSList)`
  &.MuiList-root {
    ${listCss}
    .MuiListSubheader-root {
      ${listSubheaderCss}
    }

    .MuiListItem-root {
      ${listItemCss}
      .MuiButtonBase-root {
        ${listItemButtonCss}
        ${listItemButtonSelectedCss}
        .MuiListItemIcon-root {
          ${listItemIconCss}
        }

        .MuiListItemText-root {
          ${listItemTextCss}
        }
      }
    }

    .MuiDivider-root {
      ${listItemDividerCss}
    }
  }
`;

export const NoMatches = styled(ListItem)`
  color: ${GRAY.A};
  padding: ${spacesXs}px ${spacesS}px;
`;
