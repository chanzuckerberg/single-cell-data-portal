import styled from "@emotion/styled";
import { fontBodyS, getColors } from "czifui";

export const CellCardDescription = styled.div`
  ${fontBodyS}
  font-wieght: 400;
  white-space: pre-wrap;
  ${(props) => {
    const colors = getColors(props);
    return `
      background-color: ${colors?.gray[100]};
    `;
  }}
  padding: 12px 16px 12px 16px;
  border-radius: 8px;
`;

export const Wrapper = styled.div`
  margin-top: 8px;
`;

export const Source = styled.div`
  ${fontBodyS}
  display: flex;
  justify-content: flex-end;
  margin-top: 4px;
  ${(props) => {
    const colors = getColors(props);
    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;
