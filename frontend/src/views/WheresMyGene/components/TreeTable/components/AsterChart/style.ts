import { interpolateYlOrRd } from "d3-scale-chromatic";
import styled from "styled-components";
import { SQUARE_LENGTH_PX } from "../../common/constants";

const MAX_FIRST_PART_LENGTH_PX = 16;
const MAX_SECOND_PART_LENGTH_PX = 14;

interface AsterChartProps {
  colorValue: number;
  degreeValue: number;
}

export const Container = styled.div`
  display: flex;
  align-items: center;
  justify-content: center;
  width: ${SQUARE_LENGTH_PX}px;
  height: ${SQUARE_LENGTH_PX}px;
`;

export const FirstPart = styled.div.attrs<AsterChartProps>(
  ({ colorValue, degreeValue }) => {
    const { background, height, width } = getFirstPart({
      colorValue,
      degreeValue,
    });

    return {
      style: {
        background: background,
        height: `${height}px`,
        width: `${width}px`,
      },
    };
  }
)<AsterChartProps>`
  position: relative;
  border-radius: 50%;
`;

export const SecondPart = styled.div.attrs<AsterChartProps>(
  ({ colorValue, degreeValue }) => {
    const { background, height, width, top, left } = getSecondPart({
      colorValue,
      degreeValue,
    });

    return {
      style: {
        background: background,
        height: `${height}px`,
        left: `${left}px`,
        opacity: 0.4,
        top: `${top}px`,
        width: `${width}px`,
      },
    };
  }
)<AsterChartProps>`
  position: absolute;
  border-radius: 50%;
`;

function getFirstPart({ colorValue = 0, degreeValue = 0 }) {
  const length = Math.round(MAX_FIRST_PART_LENGTH_PX * degreeValue);
  const color = interpolateYlOrRd(colorValue);
  const degree = Math.round(degreeValue * 360);

  return {
    background: `conic-gradient(${color} 0deg ${degree}deg, transparent ${degree}deg)`,
    height: length,
    width: length,
  };
}

function getSecondPart({ colorValue = 0, degreeValue = 0 }) {
  const firstPartLength = Math.round(MAX_FIRST_PART_LENGTH_PX * degreeValue);
  const length = Math.round(MAX_SECOND_PART_LENGTH_PX * degreeValue);
  const color = interpolateYlOrRd(colorValue);
  const centerOffset = firstPartLength / 2 - length / 2;

  return {
    background: color,
    height: length,
    left: centerOffset,
    top: centerOffset,
    width: length,
  };
}
