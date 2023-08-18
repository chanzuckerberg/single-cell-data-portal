import styled from "@emotion/styled";
import { spacesS } from "src/common/theme";

export const Actions = styled.div`
  align-items: center; /* required for the "more" dropdown button to center align with buttons i.e. "Download" and "Explore" */
  display: grid;
  gap: 0 ${spacesS}px;
  grid-auto-flow: column;
  justify-content: flex-start;
`;
