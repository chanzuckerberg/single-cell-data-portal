import { Button, Popover } from "@blueprintjs/core";
import { GRAY } from "src/components/common/theme";
import styled, { css } from "styled-components";
import { detailsColWidthCSS } from "../style";

export const StyledRow = styled.tr`
  align-content: center;
  vertical-align: middle;
`;

// Collection Title Column
export const CollectionTitleText = styled.a`
  font-size: 14px;
  line-height: 18px;
  letter-spacing: -0.1px;
`;
export const DOIText = styled.div`
  font-size: 12px;
  color: ${GRAY.A};
`;

// General Cell Style
export const StyledCell = styled.td`
  vertical-align: middle !important;
`;

// Details Columns
export const textClippingCSS = css`
  white-space: nowrap;
  text-overflow: clip;
  overflow: hidden;
`;
export const DetailsCell = styled(StyledCell)`
  color: ${GRAY.A};
  ${textClippingCSS}
  ${detailsColWidthCSS}
`;
export const FieldValues = styled.div`
  ${textClippingCSS}
  width: inherit;
`;
export const LeftAlignedDetailsCell = styled(DetailsCell)`
  text-align: left !important;
`;
export const RightAlignedDetailsCell = styled(DetailsCell)`
  text-align: right !important;
`;

// Field overflow popover
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
