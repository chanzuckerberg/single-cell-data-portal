import styled from "@emotion/styled";
import Image from "next/image";
import { Button } from "@czi-sds/components";
import { spacesXxxs } from "src/common/theme";

export const StyledTooltip = styled.div`
  text-align: left;
  width: 100%;

  a {
    color: inherit;
    text-decoration: underline;
  }
`;

export const TooltipButton = styled(Button)`
  font-weight: 500;
  font-size: 12px;
  min-width: 12px;
  margin-left: 2.5px;
  text-transform: none;
`;

export const StyledIconImage = styled(Image)`
  :hover {
    filter: brightness(0);
  }
`;

export const ExtraContentWrapper = styled.span`
  display: flex;
  flex-direction: row;
  align-items: center;
  gap: ${spacesXxxs}px;
  white-space: nowrap;
`;
