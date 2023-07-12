import styled from "@emotion/styled";
import { fontBodyS, getColors } from "@czi-sds/components";

export const CellGuideCardDescription = styled.div`
  ${fontBodyS}
  font-weight: 400;
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
  margin-top: 16px;
  display: flex;
  justify-content: flex-end;
  gap: 40px;
  ${(props) => {
    const colors = getColors(props);
    return `
      color: ${colors?.gray[500]};
    `;
  }}
`;

export const SourceLink = styled.div`
  white-space: nowrap;
`;
