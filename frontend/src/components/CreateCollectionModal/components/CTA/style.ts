import { GRAY, PT_GRID_SIZE_PX } from "src/components/common/theme";
import styled from "styled-components";

export const ContentWrapper = styled.div`
  color: ${GRAY.A};
  margin: ${PT_GRID_SIZE_PX}px 0 ${4 * PT_GRID_SIZE_PX}px
    ${2 * PT_GRID_SIZE_PX}px;
`;
