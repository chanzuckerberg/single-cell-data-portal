// Table.js
import styled from "@emotion/styled";
import { fontBodyS, getColors } from "czifui";

export const StyledTable = styled.table`
  width: 100%;
  table-layout: auto;
`;

export const TableWrapper = styled.div`
  max-height: 500px;
  overflow-y: auto;
`;

export const StyledHead = styled.thead`
  border-top: 0.5px solid #cccccc;
  height: 24px;
`;

export const StyledHeadCell = styled.th`
  ${fontBodyS}
  font-weight: 500;
  padding: 8px 0px 0px 0px;
  ${(props) => {
    const colors = getColors(props);
    return `color: ${colors?.gray[500]}`;
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
  padding: 4px 0px 4px 0px;
  min-width: 120px;
  word-break: break-word;
  overflow-wrap: break-word;

  a {
    display: inline-block;
    overflow-wrap: break-word;
    word-break: break-all;
  }
`;
