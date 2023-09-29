import styled from "@emotion/styled";
import { CommonThemeProps, fontBodyS, getColors } from "@czi-sds/components";
import { gray100, gray200, spacesS } from "src/common/theme";

export const StyledTable = styled.table`
  width: 100%;
  table-layout: auto;
`;

export const TableWrapper = styled.div`
  overflow-y: auto;
`;

export const StyledHead = styled.thead`
  height: 24px;
  cursor: default;
`;

export const StyledHeadCell = styled.th`
  ${fontBodyS}
  font-weight: 500;
  white-space: nowrap;
  padding-top: ${spacesS}px;
  padding-bottom: ${spacesS}px;
  ${(props) => {
    const colors = getColors(props);

    return `
    color: ${colors?.gray[500]};
    `;
  }}
`;

interface StyledRowProps extends CommonThemeProps {
  highlight: boolean;
  implementHover: boolean;
}
export const StyledRow = styled.tr<StyledRowProps>`
  background-color: ${(props) => (props.highlight ? gray100(props) : "white")};
  ${(props) =>
    props.implementHover
      ? `
    .hover-button {
      visibility: hidden;
      transition: visibility 0s ease;
    }
    &:hover .hover-button {
      visibility: visible;
    }
  `
      : ""}
  &:hover {
    background-color: ${gray200};
  }
`;

export const StyledCell = styled.td`
  ${fontBodyS}
  font-weight: 400;
  padding-top: ${spacesS}px;
  padding-bottom: ${spacesS}px;
  vertical-align: top;

  a {
    display: inline-block;
  }
`;
