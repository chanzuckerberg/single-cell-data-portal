import styled from "@emotion/styled";

export const Wrapper = styled.div`
  margin-bottom: 16px;
`;

export const Dot = styled.span`
  border-radius: 50%;
  width: 12px;
  height: 12px;

  ${({ color }: { color: string }) => {
    return `
      background-color: ${color};
    `;
  }}
`;

export const Dots = styled.div`
  display: flex;
  gap: 8px;
  margin-bottom: 10px;
`;
