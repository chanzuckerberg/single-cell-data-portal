import styled from "@emotion/styled";
import { View } from "src/views/globalStyle";
import {
  SideBarOpenButtonWrapper,
  SideBarPositioner,
} from "src/components/common/SideBar/style";
import { spacesXl, spacesXs } from "src/common/theme";

export const DatasetsView = styled(View)`
  display: grid;
  gap: ${spacesXl}px;
  place-content: flex-start stretch;
`;

export const DatasetsSideBarPositioner = styled(SideBarPositioner)`
  ${SideBarOpenButtonWrapper} {
    .MuiButton-root {
      margin-bottom: ${spacesXs}px;
    }
  }
`;
