import styled from "@emotion/styled";
import { scaleSequential } from "d3-scale";
import { interpolateMagma } from "d3-scale-chromatic";
import { maxMarkerScore } from "./constants";

const colorScale = scaleSequential(interpolateMagma)
  .domain([maxMarkerScore, 0])
  .clamp(true);

interface StyledPieProps {
  degree: number;
  size: number;
  fill: number | string;
  opacity?: number;
  center?: boolean;
}
export const StyledPie = styled.div<StyledPieProps>`
  cursor: pointer;

  ${(props) => {
    const { degree, size, fill, opacity, center } = props;
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
        opacity: ${opacity || 1};   
        position: ${center ? "relative" : "unset"};
        left: ${center ? "50%" : "unset"};
        top: ${center ? "50%" : "unset"};
        transform: ${center ? "translate(-50%, -50%)" : "unset"};
      `;
  }}
  border-radius: 50%;
`;
