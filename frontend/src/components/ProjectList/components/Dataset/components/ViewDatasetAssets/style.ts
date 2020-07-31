import { fontStyle } from "src/components/common/theme";
import { LAYOUT } from "src/components/ProjectList/common/style";
import styled from "styled-components";

export const Wrapper = styled.div`
  width: ${LAYOUT.VIEW}%;
  margin-right: ${LAYOUT.VIEW_MARGIN}%;
  display: flex;
  justify-content: center;
  align-items: center;
`;

export const SelectWrapper = styled.div`
  display: flex;
  position: relative;
`;

export const StyledSelect = styled.select`
  ${fontStyle}
  font-size: 12px;
  font-weight: 600;
  color: white;
  background-color: #318ad8;
  border: none;
  border-radius: 15px;
  appearance: none;
  padding: 0 20px;
  width: 180px;
  height: 30px;
`;

export const Arrow = styled.p`
  width: 8px;
  height: 8px;
  border: 1.5px solid #ffffff;
  border-style: none solid solid none;
  box-sizing: border-box;
  transform: rotate(45deg);
  position: absolute;
  top: calc(50% - 6px);
  right: 10%;
`;
