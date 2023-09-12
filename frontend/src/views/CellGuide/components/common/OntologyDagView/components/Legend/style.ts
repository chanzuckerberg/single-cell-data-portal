import styled from "@emotion/styled";
import {
  highlightColor,
  primaryColor,
  secondaryColor,
  tertiaryColor,
} from "../../common/constants";
import { spacesL, spacesM } from "src/common/theme";
import { Border } from "../Node/components/RectOrCircle/style";

const SQUARE_SIZE = 14;
export const LegendWrapper = styled.div`
  display: flex;
  flex-direction: row;
  column-gap: ${spacesL}px;
  padding-left: ${spacesL}px;
  padding-top: ${spacesM}px;
  position: absolute;
  z-index: 1;
  user-select: none;
  background-color: #f8f8f8;
  border-radius: 14px;
  align-items: flex-end;
`;

export const LegendItemWrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
`;

export const MarkerScoreWrapper = styled.div`
  display: flex;
  justify-content: space-between;
  width: 100%;
`;

export const LegendItem = styled.div`
  display: flex;
  justify-content: space-between;
  width: 100%;
`;

export const IsInCorpus = styled.div`
  border-radius: 50%;
  width: ${SQUARE_SIZE}px;
  height: ${SQUARE_SIZE}px;
  background-color: ${primaryColor};
`;

export const IsNotInCorpus = styled.div`
  border-radius: 50%;
  width: ${SQUARE_SIZE}px;
  height: ${SQUARE_SIZE}px;
  background-color: ${tertiaryColor};
`;

export const FlexColumn = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
`;

export const CenterText = styled.div`
  text-align: center;
`;

export const IsTargetNode = styled.div`
  border-radius: 50%;
  width: ${SQUARE_SIZE}px;
  height: ${SQUARE_SIZE}px;
  background-color: ${highlightColor};
`;

interface NoDescendantsProps {
  size: number;
}
export const NoDescendantsLine = styled.div<NoDescendantsProps>`
  position: absolute;
  width: ${(props) => props.size}px;
  height: 1px;
  background-color: ${secondaryColor};
  left: ${(props) => -props.size}px;
  top: ${(props) => props.size / 2 - 1}px;
  z-index: -1;
  border: 1px solid ${secondaryColor};
`;

export const StyledBorder = styled(Border)`
  border-left: 2px solid ${secondaryColor};
`;
