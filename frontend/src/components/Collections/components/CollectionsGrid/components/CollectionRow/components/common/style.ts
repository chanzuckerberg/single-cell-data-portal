import {
  detailsColWidthCSS,
  StyledCell,
  textClippingCSS,
} from "src/components/Collections/components/CollectionsGrid/components/common/style";
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
