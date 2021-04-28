import { Classes } from "@blueprintjs/core";
import { GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";
import { detailsColWidthCSS, textClippingCSS } from "../../../common/style";

export const StyledCell = styled.td`
  padding: ${PT_GRID_SIZE_PX * 2}px 0px !important;
  line-height: 15px;
  & > :not(a:first-child, div:first-child) {
    :not(.${Classes.TAG}) {
      display: block;
    }
    margin-top: ${PT_GRID_SIZE_PX / 2}px;
  }
`;

export const TagContainer = styled.div`
  display: flex;
  flex-direction: row;
  & > .${Classes.TAG}:first-child {
    margin-right: ${PT_GRID_SIZE_PX}px;
  }
`;

export const DetailsCell = styled(StyledCell)`
  color: ${GRAY.A} !important;
  ${textClippingCSS}
  ${detailsColWidthCSS}
  font-style: normal;
  font-weight: normal;
  font-size: 14px;
  line-height: 18px;
  letter-spacing: -0.1px;
  padding-left: ${PT_GRID_SIZE_PX * 2}px !important;
  vertical-align: middle !important;
`;

export const LeftAlignedDetailsCell = styled(DetailsCell)`
  text-align: left !important;
`;

export const StyledRow = styled.tr`
  align-content: center;
  vertical-align: top;
`;

export const RightAlignedDetailsCell = styled(DetailsCell)`
  text-align: right !important;
`;

export const ActionCell = styled(RightAlignedDetailsCell)`
  vertical-align: middle !important;
  padding-top: 0 !important;
  padding-bottom: 0 !important;
  & > .${Classes.BUTTON} {
    margin: auto 0;
  }
`;

export const ActionButtonsContainer = styled.div`
  display: flex;
  padding: 0 16px;
`;

export const ActionButton = styled.div`
  flex: 1;
`;
