import styled from "@emotion/styled";

interface StyledRectProps {
  size: number;
  fillColor: string;
}
export const StyledRect = styled.div<StyledRectProps>`
  cursor: pointer;
  height: ${(props) => props.size}px;
  width: ${(props) => props.size}px;
  background-color: ${(props) => props.fillColor};
`;
