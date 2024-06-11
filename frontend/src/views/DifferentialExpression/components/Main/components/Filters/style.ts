import styled from "@emotion/styled";
import { CommonThemeProps, getColors } from "@czi-sds/components";

export const Wrapper = styled.div<CommonThemeProps>`
  display: flex;
  flex-direction: column;
  width: 341px;
  flex-wrap: wrap;
  row-gap: 12px;
  padding: 12px;
  background-color: #f8f8f8;

  ${(props: CommonThemeProps) => {
    const colors = getColors(props);
    return `
      border: 1px solid ${colors?.gray[200]};
      border-radius: 6px;
    `;
  }}
`;

export const FlexRow = styled.div`
  display: flex;
  flex-direction: row;
  align-items: center;
`;
