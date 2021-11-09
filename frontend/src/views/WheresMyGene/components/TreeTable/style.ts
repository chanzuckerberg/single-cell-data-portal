import styled, { keyframes } from "styled-components";
import { SQUARE_LENGTH_PX } from "./common/constants";

const colorChange = keyframes`
  0%   {background-color: lightgrey; scale: 0.2;}
  25%  {background-color: darkGoldenRod; scale: 0.5;}
  50%  {background-color: darkOrange; scale: 0.8;}
  100% {background-color: red;}
`;

export const Loader = styled.div`
  width: ${SQUARE_LENGTH_PX}px;
  height: ${SQUARE_LENGTH_PX}px;
  border-radius: 50%;
  background-color: lightgrey;
  animation: ${colorChange} 3s infinite;
`;
