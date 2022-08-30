import styled from "@emotion/styled";
import { textClippingCSS } from "src/components/Collections/components/Grid/common/style";
import { GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";

export const FieldValues = styled.div`
  ${textClippingCSS}
  color: ${GRAY.A};
`;

export const ContentWrapper = styled.div`
  display: flex;
  flex-direction: row;
  padding: ${2 * PT_GRID_SIZE_PX}px;
`;

export const ContentColumn = styled.div`
  display: flex;
  flex-direction: column;
  &:not(:last-child) {
    margin-right: ${PT_GRID_SIZE_PX * 3}px;
  }
  min-width: 160px;
`;
