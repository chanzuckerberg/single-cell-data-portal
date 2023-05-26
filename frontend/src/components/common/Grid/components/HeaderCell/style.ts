import styled from "@emotion/styled";
import { ALIGNMENT } from "src/components/common/Grid/common/entities";
import { CommonThemeProps, getSpaces, getColors } from "@czi-sds/components";
import { css } from "@emotion/react";

const grey400 = (props: CommonThemeProps) => getColors(props)?.gray[400];
const primary400 = (props: CommonThemeProps) => getColors(props)?.primary[400];
const primary500 = (props: CommonThemeProps) => getColors(props)?.primary[500];
const primary600 = (props: CommonThemeProps) => getColors(props)?.primary[600];
const spacesXxs = (props: CommonThemeProps) => getSpaces(props)?.xxs;

interface Props extends CommonThemeProps {
  alignment: ALIGNMENT;
  isSortable: boolean;
  isSorted: boolean;
}

export const Header = styled("span")`
  align-items: center;
  align-self: center;
  display: flex;
  gap: ${spacesXxs}px; /* gap between header and sort icon */

  span {
    min-width: 0; /* facilitates breaking of word on columns; flex default for min width is "auto" */
  }
`;

export const SortIcon = styled.span`
  align-items: center;
  display: flex;
  justify-content: center;

  .MuiSvgIcon-root {
    color: ${grey400};
    height: 10px;
    width: 10px;
  }
`;

const sortableHeaderCss = css`
  ${Header} {
    cursor: pointer;

    &:hover {
      color: #000000;

      ${SortIcon} .MuiSvgIcon-root {
        color: inherit;
      }
    }
  }
`;

const sortedHeaderCss = (props: Props) => css`
  ${Header} {
    color: ${primary400(props)};

    &:hover {
      color: ${primary500(props)};
    }

    &:active {
      color: ${primary600(props)};
    }

    ${SortIcon} .MuiSvgIcon-root {
      color: inherit;
    }
  }
`;

export const HeaderCell = styled("th")<Props>`
  display: flex;
  gap: ${spacesXxs}px; /* gap between header and count */
  justify-content: ${(props) =>
    props.alignment === ALIGNMENT.LEFT ? "flex-start" : "flex-end"};

  ${({ isSortable }) =>
    isSortable &&
    css`
      ${sortableHeaderCss}
    `}

  ${(props) =>
    props.isSorted &&
    css`
      ${sortedHeaderCss(props)}
    `}
`;
