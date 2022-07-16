import styled, { css } from "styled-components";
import { HEADER_HEIGHT_PX } from "../Header/style";

export const Wrapper = styled.div`
  display: flex;
  flex-direction: column;
  height: 100vh; /* required for content scroll */
  min-height: 100vh; /* required for full height flex on child */
  /* TODO(seve): #2755 */
  overflow-x: hidden; /* responsive requirement; facilitates hiding of content when viewport is resized and layout min width is applied */
`;

export const contentWrapper = css`
  padding: 24px 40px;
`;

export const MainWrapper = styled.div`
  display: grid; /* required: ensures any remaining viewport height allocated to main content is observed; ancestor component heights are unspecified and so any height specification will revert to "auto" */
  flex: 1; /* sticks footer to bottom of viewport and initial render of main content is at full viewport height while data is loading. */
  margin-top: ${HEADER_HEIGHT_PX}px; /* positions content below fixed header */
`;

export const DefaultMainWrapper = styled(MainWrapper)`
  main {
    overflow: auto;
  }
`;

export const SidebarMainWrapper = styled(MainWrapper)`
  main {
    display: grid;
    grid-template-areas: "leftsidebar content rightsidebar";
    grid-template-columns: auto 1fr auto; /* grid columns for sidebar and corresponding content. */
  }
`;

export const StyledDocsLayout = styled(MainWrapper)`
  justify-content: center; /* Allows for center resizing based on viewport width */
  main {
    display: grid;
    grid-template-areas: "leftsidebar content rightsidebar";
    grid-template-columns: 368px auto 253px;
    grid-column-gap: 80px;
    max-width: 1440px;
  }
`;
