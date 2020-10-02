import { GREY } from "src/components/common/theme";
import styled from "styled-components";

export const CodeWrapper = styled.div`
  position: relative;
  width: 511px;
`;

export const Code = styled.code`
  border: 1px solid ${GREY.VERY_LIGHT};
  box-sizing: border-box;
  border-radius: 4px;
  height: 52px;
  overflow: hidden;
  display: block;
  font-size: 13px;
  line-height: 18px;
  margin-bottom: 5px;
  padding: 8px;
`;

export const CodeMask = styled.div`
  position: absolute;
  top: 0;
  width: 100%;
  height: 100%;
  opacity: 0;
  background-color: rgba(0, 118, 220, 0.9);
  color: white;
  display: flex;
  align-items: center;
  justify-content: center;
  font-size: 16px;
  cursor: pointer;

  :hover {
    opacity: 1;
  }
`;

export const Tip = styled.div`
  color: ${GREY.TIPS};
  font-size: 12px;
  width: 490px;
`;
