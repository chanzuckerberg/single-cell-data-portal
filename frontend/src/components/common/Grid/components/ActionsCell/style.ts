import styled from "styled-components";

export const Actions = styled.div`
  display: grid;
  grid-auto-columns: 32px; /* specifies an action button slot of 32px */
  grid-auto-flow: column;
  margin-top: -8px; /* visually top aligns action button icon with td padding and maintains desired action button clickable area */
`;
