import styled, { css } from "styled-components";
import { layout } from "../common/layout";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  min-height: 100vh; /* required for full height flex on child */
`;

export const contentWrapper = css`
  margin: 24px 0;
  padding: 0 40px;
`;

export const MainWrapper = styled.div`
  ${layout}
  display: grid; /* required: ensures any remaining viewport height allocated to main content is observed; ancestor component heights are unspecified and so any height specification will revert to "auto" */
  flex: 1; /* sticks footer to bottom of viewport and initial render of main content is at full viewport height while data is loading. */
`;

export const DefaultMainWrapper = styled(MainWrapper)`
  ${contentWrapper}
`;

export const SidebarMainWrapper = styled(MainWrapper)`
  main {
    display: grid;
    grid-template-areas: "sidebar content";
    grid-template-columns: auto 1fr; /* grid columns for sidebar and corresponding content. */
  }
`;
