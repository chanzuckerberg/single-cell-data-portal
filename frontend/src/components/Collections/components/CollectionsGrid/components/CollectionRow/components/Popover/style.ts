import { GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";
import { textClippingCSS } from "../../../common/style";

export const FieldValues = styled.div`
  ${textClippingCSS}
  color: ${GRAY.A};
`;

export const ContentWrapper = styled.div`
  padding: ${2 * PT_GRID_SIZE_PX}px;
`;
