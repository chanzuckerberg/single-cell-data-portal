import styled from "@emotion/styled";
import { spacesXxs, spacesS, gray400 } from "src/common/theme";

export const StyledTooltip = styled("div")`
  padding: ${spacesXxs}px 0;
  text-align: center;
  line-height: 20px;
  span {
    font-style: italic;
  }
`;

export const TooltipContent = styled("div")`
  ul,
  li {
    margin-top: 0;
    margin-bottom: 0;
    padding: 0;
    font-style: italic;
  }
`;

export const TooltipTrigger = styled("div")`
  padding-left: ${spacesS}px;
  display: inline;
  color: ${gray400};
  span {
    font-weight: 400 !important;
  }
`;
