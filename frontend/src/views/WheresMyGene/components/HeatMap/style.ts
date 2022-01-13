import styled from "styled-components";

export const Container = styled.div`
  height: 75vh;
  width: 80vw;
  overflow: scroll;
  position: relative;
`;

export const Loader = styled.div`
  position: fixed;
  top: 75px;
  left: 50vw;
  width: 200px;
  height: 50px;
  background-color: white;
  display: flex;
  align-items: center;
  gap: 10px;
  justify-content: center;
  box-shadow: 0px 0px 0px rgba(16, 22, 26, 0.1),
    0px 4px 8px rgba(16, 22, 26, 0.2);
`;
