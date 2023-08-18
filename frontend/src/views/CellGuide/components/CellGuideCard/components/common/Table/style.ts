import styled from "@emotion/styled";
import { fontBodyS, getColors, getSpaces } from "@czi-sds/components";
import { spacesM, spacesS } from "src/common/theme";

export const StyledTable = styled.table`
  width: 100%;
  table-layout: auto;
`;

export const TableWrapper = styled.div`
  max-height: 500px;
  overflow-y: auto;
`;

export const StyledHead = styled.thead`
  height: 24px;
  cursor: default;
`;

export const StyledHeadCell = styled.th`
  ${fontBodyS}
  font-weight: 500;
  ${(props) => {
    const colors = getColors(props);
    const spacings = getSpaces(props);

    return `
    color: ${colors?.gray[500]};
    padding: ${spacings?.s}px ${spacings?.m}px ${spacings?.s}px ${spacings?.m}px;`;
  }}
`;

interface StyledRowProps {
  highlight: boolean;
}
export const StyledRow = styled.tr<StyledRowProps>`
  background-color: ${(props) => (props.highlight ? "#F8F8F8" : "white")};
`;

export const StyledCell = styled.td`
  ${fontBodyS}
  font-weight: 400;
  padding: ${spacesS}px ${spacesM}px;
  min-width: 120px;
  word-break: break-word;
  overflow-wrap: break-word;
  vertical-align: top;

  a {
    display: inline-block;
    overflow-wrap: break-word;
  }
`;
