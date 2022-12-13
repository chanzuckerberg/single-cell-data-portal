import styled from "@emotion/styled";
import {
  SideBar,
  SideBarPositioner as RawSideBarPositioner,
} from "src/components/common/SideBar/style";
import { FOOTER_HEIGHT_PX } from "src/components/Footer/style";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";
import { SidebarMainWrapper } from "src/components/Layout/style";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  flex: 1;
  height: calc(100vh - ${HEADER_HEIGHT_PX}px - ${FOOTER_HEIGHT_PX}px);
`;

export const Gap = styled.div`
  height: 200px;
`;

export const Top = styled.div`
  display: flex;
  justify-content: space-between;
  margin-bottom: 20px;
`;

export const SideBarWrapper = styled(SideBar)`
  max-height: calc(100vh - ${HEADER_HEIGHT_PX}px);
`;

export const SideBarPositioner = styled(RawSideBarPositioner)`
  max-height: calc(100vh - ${HEADER_HEIGHT_PX}px);
`;

export const SideBarLayoutMainWrapper = styled(SidebarMainWrapper)`
  height: calc(100vh - ${HEADER_HEIGHT_PX}px - ${FOOTER_HEIGHT_PX}px);
  margin-top: ${HEADER_HEIGHT_PX}px;
`;
