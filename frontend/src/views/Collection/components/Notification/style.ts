import styled from "@emotion/styled";

export const StyledNotificationWrapper = styled.div`
  left: 50%;
  position: absolute;
  top: 16px;
  transform: translateX(-50%);
  z-index: 1; /* positions above Collection hero, banner */

  .MuiPaper-root {
    margin: 0;
  }
`;
