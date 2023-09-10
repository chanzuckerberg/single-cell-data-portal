import styled from "@emotion/styled";
import {
  markerGeneModeColor,
  primaryColor,
  tertiaryColor,
} from "../../common/constants";
import { spacesL, spacesM } from "src/common/theme";

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
`;

export const LegendItemWrapper = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
`;

export const LegendItem = styled.div`
  display: flex;
  justify-content: space-between;
  width: 100%;
  max-width: ${SQUARE_SIZE * 5}px;
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

export const HasDescendants = styled.div`
  border-radius: 50%;
  width: ${SQUARE_SIZE}px;
  height: ${SQUARE_SIZE}px;

  background: repeating-linear-gradient(
    -45deg,
    #999999,
    #999999 1px,
    transparent 1px,
    transparent 4px
  );
  border: 1px solid #999999;
`;

export const HasNoDescendants = styled.div`
  width: ${SQUARE_SIZE}px;
  height: ${SQUARE_SIZE}px;

  background: repeating-linear-gradient(
    -45deg,
    #999999,
    #999999 1px,
    transparent 1px,
    transparent 4px
  );
  border: 1px solid #999999;
`;

export const YesMarkergene = styled.div`
  border-radius: 50%;
  width: ${SQUARE_SIZE}px;
  height: ${SQUARE_SIZE}px;
  background-color: ${markerGeneModeColor};
`;

export const NoMarkerGene = styled.div`
  border-radius: 50%;
  width: ${SQUARE_SIZE}px;
  height: ${SQUARE_SIZE}px;
  background-color: ${tertiaryColor};
`;
