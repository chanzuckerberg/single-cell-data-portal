import styled from "@emotion/styled";
import Image from "next/image";
import { Button } from "@czi-sds/components";

export const StyledTooltip = styled.div`
  text-align: left;
  font-size: 13px;
  line-height: 20px;
  font-weight: 500;
  padding: 8px 14px;
  width: 100%;

  a {
    color: inherit;
    text-decoration: underline;
  }
`;

export const TooltipButton = styled(Button)`
  font-weight: 500;
  font-size: 12px;
  min-width: unset;
  margin-left: 2.5px;
`;

export const StyledIconImage = styled(Image)`
  :hover {
    filter: brightness(0);
  }
`;
