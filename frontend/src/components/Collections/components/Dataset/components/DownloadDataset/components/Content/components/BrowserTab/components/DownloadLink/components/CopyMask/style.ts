import styled from "@emotion/styled";
import { fontBodyS } from "@czi-sds/components";

export const CopyMask = styled.div`
  ${fontBodyS}
  align-items: center;
  background-color: rgba(0, 118, 220, 0.9);
  color: white;
  cursor: pointer;
  display: flex;
  font-size: 16px;
  height: 100%;
  justify-content: center;
  left: 0;
  letter-spacing: normal;
  line-height: inherit;
  opacity: 0;
  position: absolute;
  top: 0;
  width: 100%;

  &:hover {
    opacity: 1;
  }
`;
