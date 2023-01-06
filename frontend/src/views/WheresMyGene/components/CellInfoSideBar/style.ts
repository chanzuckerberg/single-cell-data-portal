import styled from "@emotion/styled";
import { Button, CellHeader, CellHeaderProps, getColors } from "czifui";

export const CELL_INFO_SIDEBAR_WIDTH_PX = 400;

export const ButtonContainer = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
`;

export const GeneHeaderWrapper = styled("span")`
  display: flex;
  align-items: center;
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
  margin-right: 8px;
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
