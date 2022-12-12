import { HTMLTable } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { Button } from "czifui";

export const ButtonContainer = styled.div`
  display: flex;
  flex-direction: row;
  justify-content: space-between;
`;

export const StyledHTMLTable = styled(HTMLTable)`
  & td:nth-child(3) {
    text-align: end;
  }
  & thead td {
    color: #767676 !important;
    font-weight: 500;
  }
`;

export const TooltipButton = styled(Button)`
  font-weight: 500;
  margin-right: 8px;
`;

export const CopyGenesButton = styled(Button)`
  font-weight: 500;
  font-size: 14px;
`;

export const TissueName = styled.div`
  color: #767676;
  font-weight: 500;
  font-size: 14px;
  margin-bottom: 16px;
`;
