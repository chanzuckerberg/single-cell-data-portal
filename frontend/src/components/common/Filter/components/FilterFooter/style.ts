import styled from "@emotion/styled";
import { fontBodyXs } from "@czi-sds/components";
import {
  spacesXxs,
  spacesS,
  gray500,
  gray200,
  primary400,
  primary500,
  spacesM,
} from "src/common/theme";

export const StyledTooltip = styled("div")`
  padding: ${spacesXxs}px 0;
  text-align: left;
  line-height: 20px;
  font-weight: 400;
  width: 345px;
  a {
    text-decoration: none;
    color: inherit;
    border-bottom: 1px solid;
  }
  a:hover {
    border-bottom: 1px dotted;
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

export const TooltipLink = styled("a")`
  line-height: 20px;
`;

export const TooltipTrigger = styled("span")`
  color: inherit;
  border-bottom: 1px dotted;
  transition: all 0.3s ease;
  &:hover {
    color: black;
    border-bottom: 1px solid;
  }
`;

export const FooterContent = styled("div")`
  margin: 0 ${spacesM}px 0;
  border-top: 1px solid ${gray200};
  padding: ${spacesS}px;
  ${fontBodyXs}
  color: ${gray500};
`;

export const LinkButton = styled("button")`
  background: none;
  border: none;
  color: ${primary400};
  cursor: pointer;
  font-size: inherit;
  padding: 0 0 0 1px;
  font-weight: 400;
  text-decoration: none;
  &:hover {
    color: ${primary500};
    background: none;
  }
`;
