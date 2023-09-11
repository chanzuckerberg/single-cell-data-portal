import styled from "@emotion/styled";
import { scaleSequential } from "d3-scale";
import { interpolateMagma } from "d3-scale-chromatic";

const colorScale = scaleSequential(interpolateMagma).domain([0, 4]).clamp(true);

export const StyledRect = styled.rect`
  cursor: pointer;
`;

interface StyledCircleProps {
  fillColor: number | string;
}
export const StyledCircle = styled.circle<StyledCircleProps>`
  cursor: pointer;
  fill: ${(props) =>
    typeof props.fillColor === "number"
      ? `${colorScale(props.fillColor)}80`
      : props.fillColor};
`;

interface StyledPieProps {
  degree: number;
  size: number;
  fill: number;
}
export const StyledPie = styled.div<StyledPieProps>`
  cursor: pointer;
  ${(props) => {
    const { degree, size, fill } = props;
    return `
      width: ${size}px;
      height: ${size}px;
      background: conic-gradient(
        ${colorScale(fill)} 0deg,
        ${colorScale(fill)} ${degree}deg,
        #fff0 0deg,
        #fff0 360deg
      );
    `;
  }}
  border-radius: 50%;
`;
