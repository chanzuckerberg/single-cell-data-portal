import styled from "@emotion/styled";
import { spacesXs } from "src/common/theme";
import { CommonThemeProps } from "@czi-sds/components";

interface ScrollProps extends CommonThemeProps {
  viewsToDisplay: number;
  scrollable?: boolean;
}

export const ViewsMenu = styled.span<ScrollProps>`
  display: flex;
  flex-direction: column;
  width: fit-content;
  max-width: min-content;
  flex: 1;
  padding: ${spacesXs}px;
  gap: 8px;
`;

export const ViewsMenuOptions = styled.span`
  display: flex;
  flex-direction: row;
  gap: 8px;
`;
