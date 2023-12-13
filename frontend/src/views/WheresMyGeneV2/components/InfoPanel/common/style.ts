import styled from "@emotion/styled";
import { fontBodyXxs, fontHeaderXxs } from "@czi-sds/components";
import { gray500 } from "src/common/theme";

export const LowHigh = styled.div`
  display: flex;
  justify-content: space-between;

  /* Fixes bug where numbers were wrapping in PNG export */
  white-space: nowrap;
`;

export const Header = styled.h5`
  ${fontHeaderXxs}

  margin-bottom: 10px;
`;

export const Label = styled.label`
  ${fontBodyXxs}
  font-weight: 500;
  white-space: nowrap;
  color: ${gray500};
`;

export const Content = styled.div`
  ${fontBodyXxs}

  width: 120px;
  color: ${gray500};
`;
