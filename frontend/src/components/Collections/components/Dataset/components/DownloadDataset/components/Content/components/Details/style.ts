import { fontBodyXxs, fontHeaderS } from "@czi-sds/components";
import styled from "@emotion/styled";
import { spacesS, gray100, gray500 } from "src/common/theme";

export const NoneSelected = styled.div`
  align-items: center;
  display: flex;
  flex-direction: column;
  gap: ${spacesS}px;
  justify-content: center;
  min-height: 150px;
  border-radius: 8px;
  background-color: ${gray100};
  margin-top: ${spacesS}px;
  h4 {
    ${fontHeaderS}
    font-weight: 600;
    margin-bottom: 0;
  }
  p {
    ${fontBodyXxs}
    color: ${gray500};
    margin-bottom: 0;
  }
`;
