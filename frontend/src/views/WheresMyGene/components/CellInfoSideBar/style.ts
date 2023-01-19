import styled from "@emotion/styled";
import { Button, CellHeader, CellHeaderProps, getColors } from "czifui";
import Image from "next/image";

export const CELL_INFO_SIDEBAR_WIDTH_PX = 400;

export const ButtonContainer = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
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
  font-weight: 600;
  font-size: 16px;
  line-height: 24px;
  vertical-align: middle;
`;

export const GeneCellHeader = styled(CellHeader)`
  ${(props: CellHeaderProps) => {
    const { active = false } = props;

    const colors = getColors(props);

    return `
      &:hover {
        color: ${active ? colors?.primary[500] : "black"};

        & .MuiButtonBase-root {
          /* (thuang): Maintain button hover style  */
          color: ${colors?.primary[500]}
        }
      }
    `;
  }}
`;

export const TooltipButton = styled(Button)`
  font-weight: 500;
  font-size: 16px;
  min-width: unset;
  margin: 0px 4px 0px 4px;
`;

export const CopyGenesButton = styled(Button)`
  font-weight: 500;
  font-size: 14px;
  margin-left: 10px;
`;

export const TissueName = styled.div`
  color: #767676;
  font-weight: 500;
  font-size: 14px;
  margin-bottom: 16px;
`;
