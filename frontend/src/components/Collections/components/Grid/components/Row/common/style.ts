import { GRAY } from "src/components/common/theme";
import styled from "styled-components";
import { detailsColWidthCSS, textClippingCSS } from "../../../common/style";

export const StyledCell = styled.td`
  vertical-align: middle !important;
  padding: 0 !important;
`;

export const DetailsCell = styled(StyledCell)`
  color: ${GRAY.A};
  ${textClippingCSS}
  ${detailsColWidthCSS}
`;

export const LeftAlignedDetailsCell = styled(DetailsCell)`
  text-align: left !important;
`;

export const StyledRow = styled.tr`
  align-content: center;
  vertical-align: middle;
`;

export const RightAlignedDetailsCell = styled(DetailsCell)`
  text-align: right !important;
`;
