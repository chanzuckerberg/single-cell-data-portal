import { interpolatePlasma } from "d3-scale-chromatic";
import styled from "styled-components";
import { SQUARE_LENGTH_PX } from "../../common/constants";

const MAX_FIRST_PART_LENGTH_PX = 16;
const MAX_SECOND_PART_LENGTH_PX = 14;

export const Container = styled.div`
  display: flex;
  align-items: center;
  justify-content: center;
  width: ${SQUARE_LENGTH_PX}px;
  height: ${SQUARE_LENGTH_PX}px;
`;

export const FirstPart = styled.div`
  position: relative;
  border-radius: 50%;

  ${({ proportionalExpression, relativeExpression }) => {
    const { background, height, width } = getFirstPart({
      proportionalExpression,
      relativeExpression,
    });

    return `
      background: ${background};
      height: ${height}px;
      width: ${width}px;
    `;
  }}
`;

export const SecondPart = styled.div`
  position: absolute;
  border-radius: 50%;

  ${({ proportionalExpression, relativeExpression }) => {
    const { background, height, width, top, left } = getSecondPart({
      proportionalExpression,
      relativeExpression,
    });

    return `
      background: ${background};
      height: ${height}px;
      width: ${width}px;
      top: ${top}px;
      left: ${left}px;
      opacity: 0.4;
    `;
  }}
`;

function getFirstPart({ relativeExpression = 0, proportionalExpression = 0 }) {
  const length = Math.round(MAX_FIRST_PART_LENGTH_PX * proportionalExpression);
  const color = interpolatePlasma(relativeExpression);
  const degree = Math.round(proportionalExpression * 360);

  return {
    background: `conic-gradient(${color} 0deg ${degree}deg, transparent ${degree}deg)`,
    height: length,
    width: length,
  };
}

function getSecondPart({ relativeExpression = 0, proportionalExpression = 0 }) {
  const firstPartLength = Math.round(
    MAX_FIRST_PART_LENGTH_PX * proportionalExpression
  );
  const length = Math.round(MAX_SECOND_PART_LENGTH_PX * proportionalExpression);
  const color = interpolatePlasma(relativeExpression);
  const centerOffset = firstPartLength / 2 - length / 2;

  return {
    background: color,
    height: length,
    left: centerOffset,
    top: centerOffset,
    width: length,
  };
}
