import { GRAY, LIGHT_GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";
import { detailsColWidthCSS, textClippingCSS } from "../../../common/style";

export const StyledCell = styled.td`
  padding: 0 !important;
  padding-top: ${PT_GRID_SIZE_PX * 2}px;
`;

export const DetailsCell = styled(StyledCell)`
  color: ${GRAY.A};
  ${textClippingCSS}
  ${detailsColWidthCSS}
  padding-left: ${PT_GRID_SIZE_PX * 2}px !important;
`;

export const LeftAlignedDetailsCell = styled(DetailsCell)`
  text-align: left !important;
`;

export const StyledRow = styled.tr`
  align-content: center;
  vertical-align: top;
`;

export const StyledCollectionRow = styled(StyledRow)`
  :hover {
    background-color: ${LIGHT_GRAY.E};
    cursor: pointer;
  }
`;

export const RightAlignedDetailsCell = styled(DetailsCell)`
  text-align: right !important;
`;
