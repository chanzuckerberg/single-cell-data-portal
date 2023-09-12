import { Button, fontBodyS, fontHeaderM } from "@czi-sds/components";
import styled from "@emotion/styled";
import { primary400, spacesL, spacesM } from "src/common/theme";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

export const MobileHeaderWrapper = styled.div<{ top: number }>`
  width: 100vw;
  display: flex;
  flex-direction: column;
  position: sticky;
  top: ${(props) => props.top}px;
  z-index: 2;
`;

export const MobileHeader = styled.div`
  display: flex;
  justify-content: space-between;
  height: ${HEADER_HEIGHT_PX}px;
  align-items: center;
  align-content: center;
  background-color: white;
  padding: 0 ${spacesL}px;
  gap: 16px;
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
  padding: 0px ${spacesL}px;
`;

export const SearchBarWrapper = styled.div`
  width: 100%;
`;

export const MobileSearchBarWrapper = styled.div`
  width: 100%;
  display: flex;
  flex-direction: row;
  align-items: center;
  justify-content: space-between;
  gap: ${spacesM}px;
`;

export const StyledCancelButton = styled(Button)`
  ${fontBodyS}
  font-weight: 500;
  padding: 0;
  top: 2px;
  color: ${primary400};
`;
