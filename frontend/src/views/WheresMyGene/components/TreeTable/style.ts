import styled from "styled-components";
import { SQUARE_LENGTH_PX } from "./common/constants";

export const Square = styled.div`
  width: ${SQUARE_LENGTH_PX}px;
  height: ${SQUARE_LENGTH_PX}px;
  background-color: ${(props) => props.color};
`;
