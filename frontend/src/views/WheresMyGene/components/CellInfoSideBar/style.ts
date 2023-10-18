import { HTMLTable } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { Button, fontBodyS } from "@czi-sds/components";
import Image from "next/image";
import { fontBodyXxs, getColors } from "@czi-sds/components";
import { fontWeightSemibold, gray500 } from "src/common/theme";

export const CELL_INFO_SIDEBAR_WIDTH_PX = 400;

export const ButtonContainer = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
  margin-top: 40px;
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

export const StyledMarkerGeneHeader = styled("span")`
  color: black;
  font-weight: ${fontWeightSemibold};
  font-size: 16px;
  line-height: 24px;
  vertical-align: middle;
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

export const StyledHTMLTable = styled(HTMLTable)`
  & thead td {
    color: #767676 !important;
    font-weight: 500;
    padding: 0 !important;
  }
  & td:nth-of-type(2) {
    text-align: end;
  }

  & tr:not(thead tr) {
    border-bottom: 1px solid #cccccc;
  }
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
