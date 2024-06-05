import styled from "@emotion/styled";
import { View } from "src/views/globalStyle";
import { spacesXl } from "src/common/theme";
import {
  SideBarPositioner,
  sideBarPositionerPadding,
} from "src/components/common/SideBar/style";

export const DatasetsView = styled(View)`
  display: grid;
  gap: ${spacesXl}px;
  place-content: flex-start stretch;
`;

export const DatasetsSideBarPositioner = styled(SideBarPositioner)`
  ${sideBarPositionerPadding}
`;
