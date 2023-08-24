import {
  Button,
  CommonThemeProps,
  fontBodyS,
  fontHeaderM,
  getColors,
  getSpaces,
} from "@czi-sds/components";
import styled from "@emotion/styled";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

export const MobileHeaderWrapper = styled.div`
  width: 100vw;
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
  align-content: center;
  background-color: white;
  padding: 0 16px;
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
  padding: 0px 16px;
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

  ${(props: CommonThemeProps) => {
    const spaces = getSpaces(props);

    return `
      gap: ${spaces?.m}px
    `;
  }}
`;

export const StyledCancelButton = styled(Button)`
  ${fontBodyS}
  font-weight: 500;
  padding: 0;
  top: 2px; // Centering text
  ${(props) => {
    const colors = getColors(props);

    return `
      color: ${colors?.primary[400]}
    `;
  }}
`;
