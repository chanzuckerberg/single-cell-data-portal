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

function getFirstPart({ colorValue = 0, degreeValue = 0 }) {
  const length = Math.round(MAX_FIRST_PART_LENGTH_PX * degreeValue);
  const color = interpolateYlOrRd(colorValue);

  return {
    background: color,
    height: length,
    width: length,
  };
}
