import styled from "@emotion/styled";
import { EXPANDED_WIDTH_PX } from "src/components/common/SideBar";
import { FOOTER_HEIGHT_PX } from "src/components/Footer/style";

export const Header = styled.h1`
  margin-bottom: 12px;
  font-size: 36px;
  font-weight: bold;
`;

export const Wrapper = styled.div`
  display: flex;
  position: absolute;
  top: 150px;
  z-index: 3;
  margin-left: -10px;
  width: calc(94vw - ${EXPANDED_WIDTH_PX}px);
`;

export const Content = styled.div`
  display: flex;
  justify-content: space-between;
`;

export const ColumnOne = styled.div`
  flex: 0 1 350px;

  display: flex;
  flex-direction: column;
`;

export const ColumnTwo = styled.div`
  flex: 1 2 calc(10vw);

  display: flex;
  flex-direction: column;
  justify-content: space-evenly;
`;

export const StyledStepOne = styled.div`
  flex: 1;
  margin: 10px;
  flex: 0 1 calc(76vh - ${FOOTER_HEIGHT_PX}px);
`;

export const StyledStepTwo = styled.div`
  height: 100px !important;
  margin: 10px;
`;

export const StyledStepThree = styled.div`
  flex: 1;
  margin: 10px;
`;
