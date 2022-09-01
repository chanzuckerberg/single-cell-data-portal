import styled from "@emotion/styled";

export const Actions = styled.div`
  display: grid;
  gap: 0 8px;
  grid-auto-columns: 24px; /* specifies an action button slot of 24px */
  grid-auto-flow: column;
  margin-top: -4px; /* visually top aligns action button icon with td padding and maintains desired action button clickable area */
`;
