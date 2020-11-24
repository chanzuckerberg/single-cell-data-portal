import { Button, Popover } from "@blueprintjs/core";
import styled, { css } from "styled-components";
import { textClippingCSS } from "../../../common/style";

export const FieldValues = styled.div`
  ${textClippingCSS}
  width: inherit;
`;

const popoverButtonCSS = css`
  border-radius: 3px;
  align-self: flex-end;
  padding: 4px 4px;
  margin-left: 8px;
`;

export const StyledPopover = styled(Popover)`
  ${popoverButtonCSS}
`;

export const StyledButton = styled(Button)`
  ${popoverButtonCSS}
`;
