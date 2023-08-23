import { fontHeaderM } from "@czi-sds/components";
import styled from "@emotion/styled";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

export const MobileHeaderWrapper = styled.div`
  width: 100%;
  display: flex;
  flex-direction: column;
  position: sticky;
  top: ${HEADER_HEIGHT_PX}px;
  z-index: 2;
`;

export const MobileHeader = styled.div`
  display: flex;
  justify-content: space-between;
  height: ${HEADER_HEIGHT_PX}px;
  align-items: center;
  background-color: white;
  padding: 0 8px;
  gap: 16;
`;

export const StyledTitle = styled.div`
  ${fontHeaderM}
  text-align: center;
`;

export const MobilePageNavWrapper = styled.div`
  width: 100%;
  position: absolute;
  top: ${HEADER_HEIGHT_PX}px;
  background-color: white;
  border-bottom: 1px solid lightgrey;
  padding: 0px 16px;
`;

export const MobileSearchBarWrapper = styled.div`
  width: 100%;
  padding: 0 16px;
`;
