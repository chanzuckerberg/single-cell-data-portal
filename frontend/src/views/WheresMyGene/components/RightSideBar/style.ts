import styled from "@emotion/styled";
import {
  SideBarPositioner as RawSideBarPositioner,
  SIDEBAR_BOX_SHADOW_COLOR,
} from "src/components/common/SideBar/style";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

export const RightSideBarPositioner = styled(RawSideBarPositioner)`
  max-height: ${sectionMaxheight};
  box-shadow: inset 0px 1px 0px ${SIDEBAR_BOX_SHADOW_COLOR};
`;

export const HeaderContainer = styled.div`
  display: flex;
  justify-content: space-between;
  align-items: flex-start;
`;

export const StyledTitle = styled.div`
  font-weight: 600;
  font-size: 20px;
  line-height: 24px;
`;

/**
 * Max height for the right sidebar should be 100vh - (header height).
 * Dividing the header height by two for the split right sidebar so both components
 * equal to HEADER_HEIGHT_PX so that the full length of the right sidebar does not exceed the viewport
 */
function sectionMaxheight({ maxHeight }: { maxHeight: number }) {
  return `calc(${maxHeight}vh - ${
    maxHeight === 100 ? HEADER_HEIGHT_PX : HEADER_HEIGHT_PX / 2
  }px)`;
}
