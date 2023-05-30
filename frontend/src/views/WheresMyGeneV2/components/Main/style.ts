import { Drawer } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { CommonThemeProps, getFontWeights } from "@czi-sds/components";

export const SideBarLabel = styled("span")`
  ${(props: CommonThemeProps) => {
    const fontWeights = getFontWeights(props);

    return `
      font-weight: ${fontWeights?.semibold};
    `;
  }}
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
