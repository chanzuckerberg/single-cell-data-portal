import styled from "@emotion/styled";

interface Props {
  minHeight: number;
}

export const FilterContent = styled.span<Props>`
  display: flex;
  min-height: ${({ minHeight }) => `${minHeight}px`};
`;
