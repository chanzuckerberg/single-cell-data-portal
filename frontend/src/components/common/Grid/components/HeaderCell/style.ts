import styled from "@emotion/styled";
import { fontBodyXxs, getColors } from "czifui";

export const HeaderCell = styled.span`
  display: flex;
  gap: 4px;

  span {
    min-width: 0; /* facilitates breaking of word on columns; flex default for min width is "auto" */
  }
`;

export const CountAndTotal = styled.span`
  ${fontBodyXxs}
  ${(props) => {
    const colors = getColors(props);
    return `
      background-color: ${colors?.gray[200]};
    `;
  }}
  border-radius: 20px;
  color: #000000;
  display: flex;
  flex: none;
  font-weight: 500;
  height: 20px;
  padding: 2px 8px;
`;
