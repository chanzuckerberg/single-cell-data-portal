import RightSideBar from "src/components/common/RightSideBar";
import styled from "@emotion/styled";
import { Drawer } from "@blueprintjs/core";

export const StyledRightSideBar = styled(RightSideBar)`
  /**
   * (thuang): This is to ensure the Marker Gene Right Side Bar is on top of the
   * Gene Search Bar when the viewport width is narrow enough that they overlap
   */
  z-index: 3;
`;

export const StyledSidebarDrawer = styled(Drawer)`
  .bp4-drawer-header {
    box-shadow: none;
  }
`;

export const StyledBannerContainer = styled.div`
  display: flex;
  justify-content: center;
  padding-bottom: 20px;
`;
