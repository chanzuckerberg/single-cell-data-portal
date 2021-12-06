import styled, { css } from "styled-components";
import { layout } from "../common/layout";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  height: 100vh; /* required for content scroll */
  min-height: 100vh; /* required for full height flex on child */
  overflow-x: hidden; /* responsive requirement; facilitates hiding of content when viewport is resized and layout min width is applied */
`;

export const contentWrapper = css`
  padding: 24px 40px;
`;

export const MainWrapper = styled.div`
  ${layout}
  display: grid; /* required: ensures any remaining viewport height allocated to main content is observed; ancestor component heights are unspecified and so any height specification will revert to "auto" */
  flex: 1; /* sticks footer to bottom of viewport and initial render of main content is at full viewport height while data is loading. */
  margin-top: 48px; /* positions content below fixed header */
`;

export const DefaultMainWrapper = styled(MainWrapper)`
  min-width: unset; /* overrides min width on layout */
  overflow: auto;
`;

export const SidebarMainWrapper = styled(MainWrapper)`
  min-width: unset; /* overrides min width on layout */

  main {
    display: grid;
    grid-template-areas: "sidebar content";
    grid-template-columns: auto 1fr; /* grid columns for sidebar and corresponding content. */
  }
`;

export const View = styled.div`
  ${layout}
  ${contentWrapper}
`;
