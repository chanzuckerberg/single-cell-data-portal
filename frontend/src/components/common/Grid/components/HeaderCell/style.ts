import styled from "@emotion/styled";
import { ALIGNMENT } from "src/components/common/Grid/common/entities";
import { CommonThemeProps } from "@czi-sds/components";
import { css } from "@emotion/react";
import {
  grey400,
  primary400,
  primary500,
  primary600,
  spacesXxs,
} from "src/common/theme";

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
    color: #000000;

    ${SortIcon} .MuiSvgIcon-root {
      color: ${primary400(props)};
    }

    &:hover {
      ${SortIcon} .MuiSvgIcon-root {
        color: ${primary500(props)};
      }
    }

    &:active {
      ${SortIcon} .MuiSvgIcon-root {
        color: ${primary600(props)};
      }
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
