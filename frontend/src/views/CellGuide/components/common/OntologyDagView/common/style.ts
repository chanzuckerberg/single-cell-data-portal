import styled from "@emotion/styled";
import { scaleSequential } from "d3-scale";
import { interpolateMagma } from "d3-scale-chromatic";

const colorScale = scaleSequential(interpolateMagma).domain([4, 0]).clamp(true);

interface StyledPieProps {
  degree: number;
  size: number;
  fill: number | string;
}
export const StyledPie = styled.div<StyledPieProps>`
  cursor: pointer;
  ${(props) => {
    const { degree, size, fill } = props;
    const fillColor = typeof fill === "number" ? colorScale(fill) : fill;
    return `
        width: ${size}px;
        height: ${size}px;
        background: conic-gradient(
          ${fillColor} 0deg,
          ${fillColor} ${degree}deg,
          #fff0 0deg,
          #fff0 360deg
        );
      `;
  }}
  border-radius: 50%;
`;

interface StyledCircleProps {
  fill: number | string;
  size: number;
  opacity?: number;
}
export const StyledCircle = styled.div<StyledCircleProps>`
  border-radius: 50%;
  position: relative;
  ${(props) => {
    const { size, fill, opacity } = props;
    const fillColor = typeof fill === "number" ? colorScale(fill) : fill;
    return `
        width: ${size}px;
        height: ${size}px;
        background-color: ${fillColor};
        top: 50%;
        left: 50%;
        transform: translate(-50%, -50%);    
        opacity: ${opacity || 1};    
      `;
  }}
`;
