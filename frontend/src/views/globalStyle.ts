import { contentWrapper } from "src/components/Layout/style";
import styled from "styled-components";

export const ViewGrid = styled.div`
  display: grid;
  grid-template-columns: repeat(8, 1fr);
  width: 100%;
  margin-top: 10vh;
`;

export const View = styled.div`
  ${contentWrapper}
  grid-area: content;
`;
