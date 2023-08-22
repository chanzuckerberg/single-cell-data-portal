import styled from "@emotion/styled";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

export const MobileHeaderWrapper = styled.div`
  width: 100%;
  display: flex;
  flex-direction: column;
  position: sticky;
  top: ${HEADER_HEIGHT_PX}px;
  z-index: 99999;
`;
