import styled from "@emotion/styled";

export const StyledDiv = styled.div`
  position: absolute;
  top: 0;
  bottom: 0;
  left: 0;
  right: 0;

  height: 100%;
  min-height: 100vh;
  width: 100%;
  min-width: 100vw;

  backdrop-filter: blur(10px);

  z-index: 999;
`;
