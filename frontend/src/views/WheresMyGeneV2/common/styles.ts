import { fontHeaderS } from "@czi-sds/components";
import styled from "@emotion/styled";
import { gray100, gray300, gray500 } from "src/common/theme";

interface CellProps {
  align?: boolean;
}

export const DivTable = styled.div`
  display: table;
  width: 100%;
  max-width: 600px;
  border-collapse: collapse;
`;

export const DivTableRow = styled.div`
  display: table-row;
  line-height: 24px;
  &:nth-of-type(even) {
    background-color: ${gray100};
  }
`;

export const DivTableCell = styled.div<CellProps>`
  display: table-cell;
  padding: 4px, 0px, 4px, 0px;
  text-align: ${(props) => (props.align ? "right" : "left")};
  @media (max-width: 600px) {
    display: block;
    width: 100%;
    box-sizing: border-box;

    &:not(:last-child) {
      margin-bottom: 0.625rem;
    }
  }
`;

export const DivTableHead = styled.div`
  display: table-row;
  ${fontHeaderS}
  font-weight: 500;
  color: ${gray500};
`;

export const DivTableLegend = styled.div`
  display: table-row;
  font-weight: bold;
  font-size: 1.2em;
  color: ${gray500};
  border-bottom: 1px solid ${gray300};
`;
