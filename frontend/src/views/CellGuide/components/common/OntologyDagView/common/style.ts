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
    let fillColor = typeof fill === "number" ? colorScale(fill) : fill;
    // we cannot use `opacity` directly as it does not work with `foreignObject`
    // in safari. Instead, if opacity is defined, we convert hex to rgba and set alpha
    // to opacity.
    if (opacity) {
      const hexMatch = fillColor.slice(1).match(/.{2}/g);
      const [r, g, b] = hexMatch
        ? hexMatch.map((hex) => parseInt(hex, 16))
        : [0, 0, 0];
      fillColor = `rgba(${r}, ${g}, ${b}, ${opacity})`;
    }
    return `
        width: ${size}px;
        height: ${size}px;
        background: conic-gradient(
          ${fillColor} 0deg,
          ${fillColor} ${degree}deg,
          #fff0 0deg,
          #fff0 360deg
        );
        display: ${center ? "flex" : "unset"};
        justify-content: ${center ? "center" : "unset"};
        align-items: ${center ? "center" : "unset"};
        background-color: ${center ? "#f8f8f8" : "unset"};
      `;
  }}
  border-radius: 50%;
`;
