import styled from "@emotion/styled";
import { EXPANDED_WIDTH_PX } from "src/components/common/SideBar";
import { HEADER_HEIGHT_PX } from "src/components/Header/style";

export const Header = styled.h1`
  margin-bottom: 12px;
  font-size: 36px;
  font-weight: bold;
`;

export const Wrapper = styled.div`
  display: flex;
  position: absolute;
  top: 150px;
  margin-left: -10px;
  width: calc(96vw - ${EXPANDED_WIDTH_PX}px);
  height: calc(85vh - ${HEADER_HEIGHT_PX}px);
`;

export const Content = styled.div`
  display: flex;
  justify-content: space-between;
`;

export const ColumnOne = styled.div`
  ${isHidden}

  flex: 0 1 350px;

  z-index: 10;

  display: flex;
  flex-direction: column;
`;

export const ColumnTwo = styled.div`
  flex: 1 2 calc(10vw);

  z-index: 10;

  display: flex;
  flex-direction: column;
  justify-content: space-evenly;
`;

export const StyledStepOne = styled.div`
  flex: 1 0;
  margin: 10px;
`;

export const StyledStepTwo = styled.div`
  ${isHidden}
  height: 100px !important;
  margin: 10px;
`;

export const StyledStepThree = styled.div`
  ${isHidden}
  flex: 1;
  margin: 10px;
`;

interface IsHiddenProps {
  isHidden?: boolean;
}

function isHidden({ isHidden }: IsHiddenProps) {
  if (!isHidden) return null;

  return "visibility: hidden;";
}
