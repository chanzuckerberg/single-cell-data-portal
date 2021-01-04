import {
  detailsColWidthCSS,
  StyledCell,
  textClippingCSS,
} from "src/components/Collections/components/Grid/common/style";
import { GRAY } from "src/components/common/theme";
import styled from "styled-components";

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

// Collection Title Column
export const CollectionTitleText = styled.a`
  font-size: 14px;
  line-height: 18px;
  letter-spacing: -0.1px;
`;

export const RightAlignedDetailsCell = styled(DetailsCell)`
  text-align: right !important;
`;
