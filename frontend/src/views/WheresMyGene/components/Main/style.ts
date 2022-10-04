import { Drawer } from "@blueprintjs/core";
import styled from "@emotion/styled";
import { getFontWeights } from "czifui";

export const SideBarLabel = styled("span")`
  ${(props) => {
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
