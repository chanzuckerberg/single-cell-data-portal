import styled from "@emotion/styled";
import { Button, fontBodyS, fontHeaderS } from "@czi-sds/components";
import Image from "next/image";
import { fontBodyXxs, getColors } from "@czi-sds/components";
import {
  fontWeightSemibold,
  gray500,
  gray300,
  gray100,
} from "src/common/theme";

export const CELL_INFO_SIDEBAR_WIDTH_PX = 400;

interface CellProps {
  align?: boolean;
}

export const ButtonContainer = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  margin: 40px 0px 12px 0px;
`;

export const ButtonWrapper = styled.div`
  margin: 0;
`;

export const StyledIconImage = styled(Image)`
  :hover {
    filter: brightness(0);
  }
`;

export const GeneHeaderWrapper = styled("span")`
  display: flex;
  align-items: center;
`;

export const StyledTooltip = styled("div")`
  text-align: left;
  font-size: 13px;
  line-height: 20px;
  font-weight: 500;
  padding: 8px 14px;

  a {
    text-decoration: underline;
    color: inherit;
  }
`;

export const TooltipContent = styled("div")`
  padding: 0px 0px 8px 0px;
`;

export const TooltipLink = styled("a")`
  line-height: 20px;
`;

export const StyledMarkerGeneHeader = styled("span")`
  color: black;
  font-weight: ${fontWeightSemibold};
  font-size: 16px;
  line-height: 24px;
  vertical-align: middle;
  text-transform: capitalize;
`;

export const TooltipButton = styled(Button)`
  font-weight: 500;
  font-size: 16px;
  min-width: unset;
  margin: 0px 4px 0px 4px;
`;

export const CopyGenesButton = styled(Button)`
  ${fontBodyS}
  font-weight: 500;
  margin-left: -5px;
`;

export const MarkerStrengthContainer = styled("div")`
  display: flex;
  flex-direction: row;
  gap: 10px;
  justify-content: flex-end;
`;

export const MarkerStrengthLabel = styled("div")`
  ${fontBodyXxs}
  color: ${gray500};
  line-height: unset;
`;

export const TissueName = styled.div`
  color: #767676;
  font-weight: 500;
  font-size: 14px;
  margin-bottom: 16px;
  margin-top: 4px;
`;

export const NoMarkerGenesContainer = styled("div")`
  margin-top: 16px;
  background: #f8f8f8;

  width: 100%;
  height: 120px;

  display: flex;
  flex-direction: column;

  justify-content: center;
  text-align: center;
`;

export const NoMarkerGenesHeader = styled("span")`
  ${fontBodyS}
  color: black;
  font-weight: 500;
`;

export const NoMarkerGenesDescription = styled("span")`
  ${fontBodyXxs}
  ${(props) => {
    const colors = getColors(props);
    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;

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
